#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_SLOP } from '../modules/nf-core/bedtools/slop/main'


workflow test_bcftools_mpileup_intervals {

    input = [
        [ id:'test' ], // meta map
        [ file("../test_data/tmp1.bed", checkIfExists: true), file("test_data/tmp2.bed", checkIfExists: true) ]
    ]
    fasta = file("../test_data/sizes.genome", checkIfExists: true)

    bamlist = Channel.fromPath( "../test_data/*.md.bam" )
    .map { it[-1] as String } // get only filename
    .collectFile( name: "bamlist.txt", newLine: true, sort: true )


}

workflow {
    test_bcftools_mpileup_intervals ()
}