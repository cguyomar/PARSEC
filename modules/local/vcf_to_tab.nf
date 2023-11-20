process VCF_TO_TAB {
    tag "vcf_to_bed"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'biocontainers/bedtools:2.30.0--h7d7f7ad_1' }"

    input:
    tuple val(meta), path(vcf)
   
    output:
    path  "variant_sites.bed", emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'BEGIN {OFS="\t"} !/^#/ && \$4 !~ /,/ && \$5 !~ /,/ {print \$1, \$2, \$2, \$4, \$5}' ${vcf}  \
        > variant_sites.bed
  
    """
}