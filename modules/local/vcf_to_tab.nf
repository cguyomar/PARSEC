process VCF_TO_TAB {
    tag "vcf_to_tab"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'biocontainers/bedtools:2.30.0--h7d7f7ad_1' }"

    input:
    path(vcf)
    tuple val(meta), path(intervals)
   
    output:
    path  "chunk_*.tsv", emit: positions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    awk 'BEGIN {OFS="\t"} !/^#/ && \$4 !~ /,/ && \$5 !~ /,/ {print \$1, \$2, \$4, \$5}' ${vcf}  \
        > variant_sites.tsv

    awk -F'\t' 'NR==FNR {a[\$1,\$2,\$3]=\$4; next} 
    {
        for (i in a) {
            split(i, arr, SUBSEP);

            if (\$2 >= arr[2] && \$2 <= arr[3]) {
                output_file = "chunk" a[i] ".tsv";
                print \$0 > output_file;
                break;
            }
        }
    }' genome_intervals.bed variant_sites.tsv
    """
}