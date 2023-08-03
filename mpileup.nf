#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'


workflow test_bcftools_mpileup_intervals {

    input = [
        [ id:'test' ], // meta map
        [ file("test_data/17M29921492.md.bam", checkIfExists: true), file("test_data/17M29921497.md.bam", checkIfExists: true) ],
        [file("test_data/test_interval.bed", checkIfExists: true)]
    ]
    fasta = file("test_data/sus_scrofa.fa", checkIfExists: true)
    save_mpileup = false

    BCFTOOLS_MPILEUP ( input, fasta, save_mpileup )
}

workflow {
    test_bcftools_mpileup_intervals ()
}