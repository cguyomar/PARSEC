//
// Check input samplesheet and get read channels
//

include { SAMTOOLS_VIEW_ON_INTERVAL } from '../../modules/local/samtools_view_on_interval'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { STITCH } from '../../modules/nf-core/stitch/main'
include { VCF_TO_TAB } from '../../modules/local/vcf_to_tab'

workflow IMPUTATION {
    take:
    calling_intervals       // file with all intervals together
    imputation_intervals    // channel with one interval per emission
    bams                    // [meta, [bams], bai]
    known_sites             // path (vcf)
    genome                  // [meta, fasta, fai]
    

    main:

    //Convert vcf to tab and split by chunk
    calling_intervals.collect().view()
    VCF_TO_TAB(known_sites, calling_intervals)

    VCF_TO_TAB.out.positions
    .flatten()
    .map{ it ->
        println(it)
        [
            ['id': it.simpleName],
            it
        ]
    }
    .set { positions }
    // meta, positions

    // Combine bams and intervals
    bams_with_intervals = bams
        .combine(imputation_intervals)
        .map { meta, bam, bai, meta_interval, interval -> 
            [ 
                [
                    bam_id: meta.id,
                    chunk_id: meta_interval.id,
                    id: meta.id + "_" + meta_interval.id // id specific to bam+chunk
                ], 
                bam,
                bai,
                interval
            ]
        }
    // [meta, bam, bai, interval]

    // Subsample bams on interval
    SAMTOOLS_VIEW_ON_INTERVAL(
        bams_with_intervals,
        [[],[]],
        []
    )

    // Group bams and bais
    SAMTOOLS_VIEW_ON_INTERVAL.out.bam
        .combine(SAMTOOLS_VIEW_ON_INTERVAL.out.bai, by: 0)
        .set { splitted_bams }
    // meta, bam, bai

    // Write to file a list of bamFiles per interval
    splitted_bams
        .collectFile { meta, bam, bai ->
            fname = bam[-1] as String
            [ "${meta.chunk_id}.bamList.txt", fname + '\n' ]
        }
        .map { file -> [[id:file.simpleName], file ] } // We re-add the same meta to join later
        .set { bamList }



    // Group together bams of the same chunk
    splitted_bams
        .map { meta, bam, bai -> 
            meta = [id: meta.chunk_id] // We drop bam_id
            [ meta, bam, bai ]
            }
        .groupTuple( by: [0] )
        .set { splitted_bams }


    // Join with the list of bam files per chunk
    splitted_bams
        // .groupTuple()
        // .view()
        .join(bamList)
        .set { splitted_bams } 
    // meta, [bams], [bais], bamlist

    

    // // // Stitch needs 3 inputs :
    // // // stitch
    positions.map { meta, pos -> 
        [ 
            meta,
            pos,
            [],
            [],
            pos.readLines()[0].split('\t')?.first(),  //chr
            "4", //npop
            "100"  // ngen
        ]
    }
    .set { stitch_input }


    // stitch_input = [
    //     [ id: "positions " ],
    //     positions,
    //     [],
    //     [],
    //     "1",  //chr
    //     "4", //npop
    //     "100"  // ngen
    //     ]
    // // bams
    // // def bams = [
    // //     [ id:"test_reads" ],
    // //     bams_val,
    // //     bais_val 
    // //     ]

    // // ref genome
    // def reference = [
    //     [ id:"test_reference" ],
    //     genome[0],
    //     genome[1]
    // ]

    STITCH(
        stitch_input,
        splitted_bams,
        genome,
        []
    )


    // emit:
    // bam = SAMTOOLS_VIEW_ON_INTERVAL.out.bam  // channel: [ val(meta), [ bam ] ] ??
    // versions = SAMTOOLS_VIEW_ON_INTERVAL.out. // channel: [ versions.yml ]
}
