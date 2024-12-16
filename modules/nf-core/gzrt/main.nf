process GZRT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzrt:0.8--he4a0461_0':
        'biocontainers/gzrt:0.8--he4a0461_0' }"

    input:
    tuple val(meta), path(fastqgz)

    output:
    tuple val(meta), path("*_recovered.fastq.gz"), emit: fastqrecovered
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (fastqgz.extension != "gz") {
        error "GZRT works with .gz files only."
    }

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzrecover -o ${prefix}_recovered.fastq ${fastqgz}
    gzip ${prefix}_recovered.fastq

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"
    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_recovered.fastq.gz

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"

    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """
}
