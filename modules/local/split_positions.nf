process SPLIT_POSITIONS {
    tag "vcf_to_tab"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'biocontainers/bedtools:2.30.0--h7d7f7ad_1' }"

    input:
    tuple val(meta), val(interval), path(bed)
   
    output:
    tuple val(meta), path("res/chunk_*.tsv"), emit: positions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir res
    echo '${interval}' > int.bed
   
    chunk_id=\$(head -1 int.bed | cut -f4)
    echo chunk\$chunk_id
    bedtools intersect -a ${bed} -b int.bed | awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$4,\$5}' > res/chunk_\$chunk_id.tsv

    """
}