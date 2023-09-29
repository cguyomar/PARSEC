//
// Split bams and run Stitch
//

include { SAMTOOLS_VIEW_ON_INTERVAL } from '../../modules/local/samtools_view_on_interval'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { STITCH } from '../../modules/nf-core/stitch/main'
include { VCF_TO_TAB } from '../../modules/local/vcf_to_tab'
include { SPLIT_POSITIONS } from '../../modules/local/split_positions'
include { SAMTOOLS_SPLIT_BAM } from "../../modules/local/samtools_split_bam"

workflow IMPUTATION {
    take:
    calling_intervals       // file with all intervals together
    imputation_intervals    // file with all intervals together
    bams                    // [meta, [bams], [bais]]
    known_sites             // path (vcf)
    genome                  // [meta, fasta, fai]
    

    main:


    SAMTOOLS_SPLIT_BAM(
        bams,
        imputation_intervals,
        [[],[]]
    )

    SAMTOOLS_SPLIT_BAM.out.indexed_bams
        .flatMap { meta, bams, bais -> 
            res = []
            bams.eachWithIndex { bam, index ->
                filename = bam.getName()
                chunk_id = "chunk" + bam.simpleName
                bam_name = filename.split("\\.")
                bam_name = bam_name[1..-2].join(".")
                // meta.bam_id = meta.id
                // meta.chunk_id = chunk_id
                // meta.id = bam_name + ":" + chunk_id
                meta = [ id: chunk_id ]
                res.add([ meta, bam, bais[index] ] )
            }
            res
        }
        .set { indexed_bams_per_interval }
    // meta, bam, bai  

    // Write to file a list of bamFiles per interval
    indexed_bams_per_interval
        .collectFile { meta, bam, bai ->
            fname = bam[-1] as String
            [ "${meta.id}.bamList.txt", fname + '\n' ]
        }
        .map { file -> [[id:file.simpleName], file ] } // We re-add the same meta to join later
        .set { bamList }

    indexed_bams_per_interval
        .groupTuple( by: [0] )
        .set { indexed_bams_per_interval }

    // Join with the list of bam files per chunk
    indexed_bams_per_interval
        .join(bamList)
        .set { indexed_bams_per_interval } 
    // meta, [bams], [bais], bamlist



    // What we want : meta, [bams], [bais], bamlist

    //Convert vcf to tab and split by chunk

    // calling_intervals.view()
    VCF_TO_TAB( known_sites )

    calling_intervals.combine( VCF_TO_TAB.out.bed )
        .set { sites_to_split }
    SPLIT_POSITIONS( sites_to_split )

    SPLIT_POSITIONS.out.positions
        .map{ meta, it ->
            [
                ['id': it.simpleName],
                it
            ]
        }
        .set { positions }
        // meta, positions

    positions.filter { meta, positions -> 
        positions.countLines() > 0
    } .set { positions }

    // Prepare stitch input
    positions.map { meta, pos -> 
        [ 
            meta,
            pos,
            [],
            [],
            pos.readLines()[0].split('\t')?.first(),  //chr
            params.npop, //npop
            params.ngen  // ngen
        ]
    }
    .set { stitch_input }

    stitch_input
        .join(indexed_bams_per_interval)
        .multiMap { meta, pos,  arg1, arg2, chr, npop, ngen, bam, bai, bamlist -> 
            input1: [meta, pos, arg1, arg2, chr, npop, ngen]
            input2: [meta, bam, bai, bamlist]
        }
        .set { stitch_input }

    STITCH(
        stitch_input.input1,
        stitch_input.input2,
        genome,
        []
    )


    // emit:
    // bam = SAMTOOLS_VIEW_ON_INTERVAL.out.bam  // channel: [ val(meta), [ bam ] ] ??
    // versions = SAMTOOLS_VIEW_ON_INTERVAL.out. // channel: [ versions.yml ]
}
