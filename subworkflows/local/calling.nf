//
// Split bams according to intervals and run mpileup
//

include { SAMTOOLS_MERGE_ON_INTERVAL } from '../../modules/local/samtools_merge_on_interval'
include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main'

workflow CALLING {
    take:
    calling_intervals    // channel with one interval per emission
    bams                    // [meta, [bams], bai]
   

    main:
    // Combine bams and intervals
    bams_with_intervals = bams
        .combine(calling_intervals)
        .map { meta, bam, bai, meta_interval, interval -> 
            [ 
                meta_interval,  // bam metadata are dropped
                bam,
                bai,
                interval
            ]
        }
    // [meta, bam, bai, interval]

    // Group bams of the same interval
    bam_groupped_by_interval = bams_with_intervals
    .view()
        .groupTuple(by: [0,3])
        .map { meta, bamlist, bailist, interval ->
            [meta, bamlist, bailist, interval]
        }
    // [meta, [bams], [bais], interval]

    ///
    /// MODULE: Run Samtools merge
    ///
    SAMTOOLS_MERGE_ON_INTERVAL(
        bam_groupped_by_interval,    
        [[],[]],
        [[],[]]
        )

    ///
    /// MODULE: Run Mpileup
    ///
    BCFTOOLS_MPILEUP(
        SAMTOOLS_MERGE_ON_INTERVAL.out.bam,
        params.fasta,
        []
        )


    vcf_by_interval = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map { it ->
            [
                [id: "all_samples"],
                it[1],
                it[2]
            ]
        }
        .groupTuple()

    ///
    /// MODULE: Run Bcftools concat
    ///
    vcf_by_interval.view()
    BCFTOOLS_CONCAT(vcf_by_interval)


    ///
    /// MODULE: Run Bcftools sort
    ///
    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)

    emit:
    vcf = BCFTOOLS_SORT.out.vcf  
    // versions = SAMTOOLS_VIEW_ON_INTERVAL.out. // channel: [ versions.yml ]
}
