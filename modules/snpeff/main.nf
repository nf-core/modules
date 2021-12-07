process SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::snpeff=5.0" : null)
    container "${ task.ext.container_tag ?
            "nfcore/snpeff:${task.ext.container_tag}" :
            workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1' :
            'quay.io/biocontainers/snpeff:5.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    val   db
    path  cache

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "*.csv"                      , emit: report
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    """
    snpEff \\
        -Xmx${avail_mem}g \\
        $db \\
        $args \\
        -csvStats ${prefix}.csv \\
        $cache_command \\
        $vcf \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
