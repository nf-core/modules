process JELLYFISH_DUMP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmer-jellyfish:2.3.1--py310h184ae93_5':
        'biocontainers/kmer-jellyfish:2.3.1--py310h184ae93_5' }"

    input:
    tuple val(meta), path(jf)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: output
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("-c") ? "txt" : "fa"
    """
    jellyfish \\
        dump \\
        $args \\
        ${jf} > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version |& sed '1!d ; s/jellyfish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("-c") ? "txt" : "fa"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version |& sed '1!d ; s/jellyfish //')
    END_VERSIONS
    """
}
