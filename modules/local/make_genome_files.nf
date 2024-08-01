process MAKE_GENOME_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"
    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.chrom.sizes"), emit: chrom_sizes
    tuple val(meta), path("*.chrom.bed"), emit: chrom_bed

    script:
    """
    awk -v OFS='\t' {'print \$1,\$2'} $fai > ${meta.id}.chrom.sizes
    awk -v OFS='\t' {'print \$1,0,\$2,\$1'} $fai > ${meta.id}.chrom.bed

    """
}