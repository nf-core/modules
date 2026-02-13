process JELLYFISH_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmer-jellyfish:2.3.1--py310h184ae93_5':
        'biocontainers/kmer-jellyfish:2.3.1--py310h184ae93_5' }"

    input:
    tuple val(meta), path(fasta)
    val kmer_length
    val size

    output:
    tuple val(meta), path("${prefix}.jf"), emit: jf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    jellyfish \\
        count \\
        $args \\
        -m ${kmer_length} \\
        -s ${size} \\
        -t $task.cpus \\
        -o ${prefix}.jf \\
        ${fasta}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version |& sed '1!d ; s/jellyfish //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.jf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jellyfish: \$(jellyfish --version |& sed '1!d ; s/jellyfish //')
    END_VERSIONS
    """
}
