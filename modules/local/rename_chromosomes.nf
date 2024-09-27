process RENAME_CHROMOSOMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"
    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path(bed), emit: bed

    script:
    """
    sed -i -e 's/\\.//g' $bed 
    
    """
}