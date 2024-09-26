//
// Split bams and run Stitch
//

include { SAMTOOLS_VIEW_ON_INTERVAL } from '../../modules/local/samtools_view_on_interval'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { STITCH } from '../../modules/nf-core/stitch/main'
include { GLIMPSE_PHASE } from '../../modules/nf-core/glimpse/phase/main'
include { GLIMPSE_LIGATE } from '../../modules/nf-core/glimpse/ligate/main'
include { BEAGLE4_BEAGLE } from '../../modules/local/beagle4'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main' 
include { BCFTOOLS_VIEW } from '../../modules/nf-core/bcftools/view/main'                                   
include { VCF_TO_TAB } from '../../modules/local/vcf_to_tab'
include { SPLIT_POSITIONS } from '../../modules/local/split_positions'
include { SAMTOOLS_SPLIT_BAM } from "../../modules/local/samtools_split_bam"
include { TABIX_TABIX as TABIX_STITCH } from '../../modules/nf-core/tabix/tabix/main'   
include { TABIX_TABIX as TABIX_BEAGLE } from '../../modules/nf-core/tabix/tabix/main'   
include { TABIX_TABIX as TABIX_GLIMPSE } from '../../modules/nf-core/tabix/tabix/main'   
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
                
                // Regexp to capture chunk name ("chr_i") in file name
                matcher = (bam.name =~ /^(([^._]+_)?[^.]+)(?:\..*)?/)
                chunk_id = "chunk_" +  matcher[0][1]
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

    print(indexed_bams_per_interval)

    // What we want : meta, [bams], [bais], bamlist

    //Convert vcf to tab and split by chunk

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
                [], // input
                [], // rdata
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
            genome.collect(),
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

    intervals = imputation_intervals.join(calling_intervals)

    // Turn bed like intervals in chr:start-end
    intervals.map  {meta, int1, int2 -> 
        split1 = int1.split("\t")
        res1 = "${split1[0]}:${split1[1]}-${split1[2]}"
        split2 = int2.split("\t")
        res2 = "${split2[0]}:${split2[1]}-${split2[2]}"
        return([meta, res1, res2])
    }.set { intervals }  // meta, calling_interval, imputation_interval
    
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
        
        ///
        /// Glimpse Phase
        /// needs val(meta) , path(input), path(input_index), path(samples_file), val(input_region), val(output_region), path(reference), path(reference_index), path(map)
        GLIMPSE_PHASE(glimpse_input)


        ///
        /// Index Glimpse output
        ///
        TABIX_GLIMPSE( GLIMPSE_PHASE.out.phased_variant )

        GLIMPSE_PHASE.out.phased_variant
            .join(TABIX_GLIMPSE.out.tbi)
            .set { indexed_glimpse_output }

        ligate_input = indexed_glimpse_output
            .map { it -> [[id:"glimpse_output"], it[1], it[2]]}
            .groupTuple( by: 0 )
        ///
        /// Glimpse Ligate
        ///
        GLIMPSE_LIGATE( ligate_input )
    }


    if (params.imputation_tool == "beagle4"){

        // Beagle is run on the calling intervals, and vcf is then subsetted on the imputation intervals
        intervals.map { meta, calling_int, imputation_int  -> 
            [
                meta,
                calling_int
            ]}.set{ beagle_intervals }
        beagle_intervals.view()

        if (ref_panel==null){
            BEAGLE4_BEAGLE(
            sparse_vcf.collect(), // allows to run several times with only one vcf
            beagle_intervals,
            [],
            [],
            [],
            [],
        )
        } else {
            BEAGLE4_BEAGLE(
            sparse_vcf.collect(), // allows to run several times with only one vcf
            beagle_intervals,
            ref_panel,
            [],
            [],
            [],
        )
        }


        

        // bcftools index
        BCFTOOLS_INDEX(BEAGLE4_BEAGLE.out.vcf)

        BEAGLE4_BEAGLE.out.vcf
            .join(BCFTOOLS_INDEX.out.csi, by: 0)
            .set { indexed_vcfs }
        // meta, vcf, tbi

        // Write intervals to files as the bcftools view module does not accept vals
        calling_intervals.collectFile(sort: false) { it -> 
            [ it[0].id + ".bed", it[1] ]
        }.merge(calling_intervals) // bed, meta, interval
        .map { it -> [ it[1], it[0] ] }
        .set { calling_intervals_as_bed } // meta, bed

        // bcftools view
        indexed_vcfs
        .join(calling_intervals_as_bed)
            .multiMap { it -> 
                intervals: it[3]
                vcf: [ it[0], it[1], it[2] ]
            }
            .set { indexed_vcf_for_view }  

        BCFTOOLS_VIEW(
            indexed_vcf_for_view.vcf,
            indexed_vcf_for_view.intervals,
            [],
            []
        )


        TABIX_BEAGLE( BCFTOOLS_VIEW.out.vcf )

        TABIX_BEAGLE.out.tbi
        .join(BEAGLE4_BEAGLE.out.vcf)
        .map{ it -> [ [id: "imputed"], it[2], it[1] ] }
        .groupTuple() 
        .set { beagle_vcf_indexed }
        // meta, [vcf.gz], [tbi]
        
         
        // MODULE : BCFTOOLS CONCAT
        BCFTOOLS_CONCAT( beagle_vcf_indexed.collect() )
    
    }

    // emit:
    // bam = BCFTOOLS_CONCAT.out.vcf  // channel: [ val(meta), vcf ] ??
    // versions = SAMTOOLS_VIEW_ON_INTERVAL.out. // channel: [ versions.yml ]
}
