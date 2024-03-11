process BEAGLE4_BEAGLE {
    tag "$meta_vcf.id" + ":" + "$meta_interval.id"
    label 'process_high'

    conda "bioconda::beagle=4.1_21Jan17.6cc.jar"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/beagle:4.1_21Jan17.6cc.jar--0 ' :
        'biocontainers/beagle:4.1_21Jan17.6cc.jar--0 ' }"

    input:
    tuple val(meta_vcf), path(vcf)
    tuple val(meta_interval), val(interval)
    path(refpanel)
    path(genmap)
    path(exclsamples)
    path(exclmarkers)

    output:
    tuple val(meta_interval), path("*.vcf.gz")     , emit: vcf
    tuple val(meta_interval), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_interval.id}.bglout"
    def ref_command = refpanel ? "ref=$refpanel" : ""
    def map_command = genmap ? "map=$genmap" : ""
    def interval_command = interval ? "chrom=$interval" : ""
    def excludesamples_command = exclsamples ? "excludesamples=$exclsamples" : ""
    def excludemarkers_command = exclmarkers ? "excludemarkers=$exclmarkers" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[beagle] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    beagle -Xmx${avail_mem}M \\
        gl=${vcf} \\
        out=${prefix} \\
        $args \\
        ${interval_command} \\
        ${ref_command} \\
        ${map_command} \\
        ${excludesamples_command} \\
        ${excludemarkers_command} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: \$(beagle 2>&1 |head -n1 | sed -rn 's/beagle\\.(.*)\\.jar \\(version (.*)\\)/\\2rev\\1/p')
    END_VERSIONS
    """
}