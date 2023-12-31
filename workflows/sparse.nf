/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSparse.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { IMPUTATION } from '../subworkflows/local/imputation'
include { CALLING } from '../subworkflows/local/calling'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { QUALIMAP_BAMQC } from '../modules/nf-core/qualimap/bamqc/main'  
include { BEDTOOLS_MAKEWINDOWS } from '../modules/nf-core/bedtools/makewindows/main' 
include { BEDTOOLS_SLOP } from '../modules/nf-core/bedtools/slop/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/// 
/// MODULE: Local modules
///
include { SAMTOOLS_MERGE_ON_INTERVAL } from '../modules/local/samtools_merge_on_interval'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SPARSE {

    ch_versions = Channel.empty()


    // Validate input parameters
    if (params.imputation_tool != "stitch" & params.imputation_tool != "glimpse") {
        exit 1, 'Imputation tool should be one of stich, glimpse'
    }

    if (params.imputation_tool == "glimpse" & !params.ref_panel) {
        exit 1, 'Imputation tool Glimpse requires a reference panel supplied with --ref_pannel'
    } else if (params.imputation_tool=="glimpse") {
        reference_panel = 
            [
                [ id:"reference_panel" ],
                file(params.ref_panel, checkIfExists: true),
            ]
        
    } else {
        reference_panel = null
    }

    // //
    // // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    // //
    // INPUT_CHECK (
    //     file(params.input)
    // )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema
    bam_channel = Channel.fromPath(
        params.bam,
        checkIfExists: true
        ).map {
        it ->
        [
            ['id': it.baseName], // meta map
            it
        ]
    }

    ///
    /// MODULE: Run Samtools index
    ///
    SAMTOOLS_INDEX(
        bam_channel
    )

    bam_channel.combine(SAMTOOLS_INDEX.out.bai, by: 0)
        .set { indexed_bams }

    ///
    /// MODULE: Run MakeWindows
    ///

    BEDTOOLS_MAKEWINDOWS(
        [
            [ id:'genome_intervals'],
            file(params.genome_bed, checkIfExists: true)
        ]
    )

    BEDTOOLS_SLOP(
        BEDTOOLS_MAKEWINDOWS.out.bed,
        file(params.genome_sizes, checkIfExists: true)
    )

    BEDTOOLS_MAKEWINDOWS.out.bed
        .join(BEDTOOLS_SLOP.out.bed)
        .set { bedFiles }

    // Turn interval files into one bed file per line/chunk
    bedFiles
        .flatMap { meta,bed,sloppedBed ->
            lines = bed.readLines()
            lines_slopped = sloppedBed.readLines()
            chunk_ids = []
            res = []
            def nb_interval = 1
            for (int i = 0; i < lines.size(); i++) {
                def line = lines[i]
                def line_slopped = lines_slopped[i]
                chunk_id = "chunk_" + nb_interval.toString()
                chunk_ids.add(chunk_id)
                nb_interval += 1

                res.add([chunk_id,line,line_slopped])
            }
            res
        }
        .multiMap { it ->
            calling: [ it[0], it[1] ]
            imputation: [ it[0], it[2] ]
        }
        .set { intervals }
        
    
    /// Read interval ids to update the meta map
    /// [interval, chunk_id]
    intervals.calling
        .map {
            it -> [[ id:it[0]],it[1]]
        }
        .set { intervals_for_calling }

    intervals.imputation
        .map {
            it -> [[ id:it[0]],it[1]]
        }
        .set { intervals_for_imputation }


       if ( !params.sparse_variants )  {
        CALLING(
            intervals_for_calling_as_bed,
            indexed_bams
        )
        CALLING.out.vcf
            .set { sparse_variants }
    } else {
        sparse_variants = [
            [ id:"sparse_variants" ],
            file(params.sparse_variants)
        ]
    }

    def reference_genome = [
        [ id:"test_reference" ],
        file(params.fasta, checkIfExists: true),
        file(params.fai, checkIfExists: true),
    ]

    // reference_genome[0].view()
  
    IMPUTATION(
        intervals_for_calling,
        intervals_for_imputation,
        indexed_bams,
        Channel.fromList([sparse_variants]),
        Channel.fromList([reference_genome]),
        Channel.fromList([reference_panel])
    )
    // //
    // // MODULE: Run BamQC
    // //
    // QUALIMAP_BAMQC (
    //     bam_channel,
    //     []
    // )
    // ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

//     CUSTOM_DUMPSOFTWAREVERSIONS (
//         ch_versions.unique().collectFile(name: 'collated_versions.yml')
//     )

//     //
//     // MODULE: MultiQC
//     //
//     workflow_summary    = WorkflowSparse.paramsSummaryMultiqc(workflow, summary_params)
//     ch_workflow_summary = Channel.value(workflow_summary)

//     methods_description    = WorkflowSparse.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
//     ch_methods_description = Channel.value(methods_description)

//     ch_multiqc_files = Channel.empty()
//     ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//     // ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.zip.collect{it[1]}.ifEmpty([]))

//     MULTIQC (
//         ch_multiqc_files.collect(),
//         ch_multiqc_config.toList(),
//         ch_multiqc_custom_config.toList(),
//         ch_multiqc_logo.toList()
//     )
//     multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
