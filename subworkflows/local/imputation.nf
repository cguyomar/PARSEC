//
// Split bams and run Stitch
//

include { SAMTOOLS_VIEW_ON_INTERVAL } from '../../modules/local/samtools_view_on_interval'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { STITCH } from '../../modules/nf-core/stitch/main'
include { GLIMPSE_PHASE } from '../../modules/nf-core/glimpse/phase/main'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { VCF_TO_TAB } from '../../modules/local/vcf_to_tab'
include { SPLIT_POSITIONS } from '../../modules/local/split_positions'
include { SAMTOOLS_SPLIT_BAM } from "../../modules/local/samtools_split_bam"
include { TABIX_TABIX as TABIX_STITCH } from '../../modules/nf-core/tabix/tabix/main'   
include { TABIX_TABIX as TABIX_REF } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_SPARSE } from '../../modules/nf-core/tabix/tabix/main'   

      

workflow IMPUTATION {
    take:
    calling_intervals       // chanel emitting each interval as a [meta, val]
    imputation_intervals    // chanel emitting each interval as a [meta, val]
    bams                    // [meta, [bams], [bais]]
    sparse_vcf             // [meta, path (vcf)]
    genome                  // [meta, fasta, fai]
    ref_panel               // [meta, path(vcf)]
    

    main:


    calling_intervals.map { it ->
            it[1]
        }
        .collectFile(name: "calling_intervals.bed", newLine: true)
        .map { it -> [
            [ id:'calling_intervals' ],
            it
            ]
        }
        .set{ calling_intervals_bed }

    imputation_intervals.map { it ->
            it[1]
        }
        .collectFile(name: "imputation_intervals.bed", newLine: true)
        .map { it -> [
            [ id:'imputation_intervals' ],
            it
            ]
        }
        .set{ imputation_intervals_bed }
    bams.combine(imputation_intervals_bed)
        .set{ bams_with_intervals }

    SAMTOOLS_SPLIT_BAM(
        bams_with_intervals,
        [[],[]]
    )

    SAMTOOLS_SPLIT_BAM.out.indexed_bams
        .flatMap { meta, bams, bais -> 
            res = []
            bams.eachWithIndex { bam, index ->
                filename = bam.getName()
                chunk_id = "chunk_" + bam.simpleName
                // bam_name = filename.split("\\.")
                // bam_name = bam_name[1..-2].join(".")
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
        .collectFile(sort: true) { meta, bam, bai ->
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
    VCF_TO_TAB( sparse_vcf )

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
        positions.countLines() > 10
    } .set { positions }

    if (params.imputation_tool == "stitch"){
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

    // MODULE : BGZIP TABIX
        TABIX_STITCH( STITCH.out.vcf )

        TABIX_STITCH.out.tbi
        .join(STITCH.out.vcf)
        .map{ it -> [ [id: "imputed"], it[2], it[1] ] }
        .groupTuple() 
        .set { stitch_vcf_indexed }
        // meta, [vcf.gz], [tbi]
        
         
    // MODULE : BCFTOOLS CONCAT
    BCFTOOLS_CONCAT( stitch_vcf_indexed.collect() )
    }
    
    if (params.imputation_tool == "glimpse"){
         // Prepare data for Glimpse

        // Index sparse vcf
        TABIX_SPARSE( sparse_vcf )

        sparse_vcf
            .join( TABIX_SPARSE.out.tbi )
            .set { indexed_sparse_vcf }

        // Index reference vcf
        TABIX_REF( ref_panel )

        ref_panel
            .join( TABIX_REF.out.tbi )
            .set { indexed_ref_panel }

        intervals = imputation_intervals.join(calling_intervals)

        // Turn bed like intervals in chr:start-end
        intervals.map  {meta, int1, int2 -> 
            split1 = int1.split("\t")
            res1 = "${split1[0]}:${split1[1]}-${split1[2]}"
            split2 = int2.split("\t")
            res2 = "${split2[0]}:${split2[1]}-${split2[2]}"
            return([meta, res1, res2])
        }.set { intervals }
        
        empty_ch = Channel.of([[]])
        glimpse_input = indexed_sparse_vcf
            .combine(empty_ch) // samples file
            .combine(intervals)
            .combine(indexed_ref_panel)
            .combine(empty_ch) // map
            .map { meta_vcf, input_vcf, input_vcf_idx, samples_file, meta_intervals, input_interval, output_interval, meta_reference, reference, reference_idx, map ->
                [
                    meta_intervals,
                    input_vcf,
                    input_vcf_idx,
                    samples_file,
                    input_interval,
                    output_interval,
                    reference,
                    reference_idx,
                    map
                ]
            }


        // needs val(meta) , path(input), path(input_index), path(samples_file), val(input_region), val(output_region), path(reference), path(reference_index), path(map)
        GLIMPSE_PHASE(glimpse_input)
    }

    // emit:
    // bam = BCFTOOLS_CONCAT.out.vcf  // channel: [ val(meta), vcf ] ??
    // versions = SAMTOOLS_VIEW_ON_INTERVAL.out. // channel: [ versions.yml ]
}
