process SAMTOOLS_SPLIT_BAM {
    tag "$meta_bam.id" + ":" + "$meta_intervals.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta_bam), path(input), path(index), val(meta_intervals), path(intervals)
    tuple val(meta_fasta), path(fasta)

    output:
    tuple val(meta_bam), path("bam_split/*.bam"), path("bam_split/*.bai") , optional:true, emit: indexed_bams
    // tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    // tuple val(meta), path("*.csi")         , optional:true, emit: csi
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta_bam.id}"
    // def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    mkdir -p bam_split

    while IFS=\$'\t' read -r chrom start end name
    do
        #output_bam="bam_split/\$chrom_\$start_\$end.bam"
        samtools view -b "${input}" "\$chrom:\$start-\$end" > bam_split/\$name.${input}
        samtools index bam_split/\$name.${input}
    done < ${intervals}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
