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
    path  "res/chunk_*.tsv", emit: positions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir res
    awk 'BEGIN {OFS="\t"} !/^#/ && \$4 !~ /,/ && \$5 !~ /,/ {print \$1, \$2, \$2, \$4, \$5}' ${vcf}  \
        > variant_sites.bed

    # split bed to make thinks faster
    for chr in `cut -f 1 ${intervals} | sort | uniq`;
    do
        echo \$chr
        grep -w \$chr variant_sites.bed > \$chr.bed
    done

    while read line 
    do 
        chunk_id=\$(echo \$line | cut -d " " -f4)
        chr=\$(echo \$line | cut -d " " -f1)
        echo chunk\$chunk_id
        echo \$line | sed -e 's/ /\t/g' > int.bed
        bedtools intersect -sorted -a \$chr.bed -b int.bed | awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$4,\$5}' > res/chunk\$chunk_id.tsv
    done < ${intervals}

    """
}