process GZRT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzrt:0.9.1--h577a1d6_0':
        'biocontainers/gzrt:0.9.1--h577a1d6_0' }"

    input:
    tuple val(meta), path(fastqgz)

    output:
    tuple val(meta), path("${prefix}.fastq.gz"), emit: recovered
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (fastqgz.extension != "gz") {
        error "GZRT works with .gz files only."
    }

    prefix = task.ext.prefix ?: "${meta.id}_recovered"
    if ("${fastqgz}" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    gzrecover -p ${fastqgz} | gzip > ${prefix}.fastq.gz

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"
    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_recovered"
    if ("${fastqgz}" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    echo "" | gzip > ${prefix}.fastq.gz

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"

    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """
}
