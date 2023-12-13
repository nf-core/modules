
process TRUVARI_CONSISTENCY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:4.1.0--pyhdfd78af_0':
        'biocontainers/truvari:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf1), path(vcf2), path(vcf3)

    output:
    tuple val(meta), path("*.consistency") , emit: consistency
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf3   = vcf3 ?: "$vcf3"

    """
    truvari \\
        consistency \\
        $args \\
        $vcf1 \\
        $vcf2  \\
        $vcf3 > ${prefix}.consistency

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consistency

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' ))
    END_VERSIONS
    """
}
