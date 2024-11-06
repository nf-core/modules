process PBSV_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbsv:2.9.0--h9ee0642_0':
        'biocontainers/pbsv:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta),  path(svsig)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbsv \\
        call \\
        $args \\
        -j ${task.cpus} \\
        ${fasta} \\
        ${svsig} \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsv: \$(pbsv --version |& sed '1!d ; s/pbsv //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsv: \$(pbsv --version |& sed '1!d ; s/pbsv //')
    END_VERSIONS
    """
}
