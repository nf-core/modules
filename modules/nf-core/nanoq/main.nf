process NANOQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2'
        : 'biocontainers/nanoq:0.10.0--h031d066_2'}"

    input:
    tuple val(meta), path(ontreads)
    val(output_format) //One of the following: fastq, fastq.gz, fastq.bz2, fastq.lzma, fasta, fasta.gz, fasta.bz2, fasta.lzma.

    output:
    tuple val(meta), path("*.{stats,json}")            , emit: stats
    tuple val(meta), path("${prefix}.${output_format}"), emit: reads
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    nanoq -i ${ontreads} \\
        ${args} \\
        -r ${prefix}.stats \\
        -o ${prefix}.${output_format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    echo "" | gzip > ${prefix}.${output_format}
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """
}
