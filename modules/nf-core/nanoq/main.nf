process NANOQ {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2' :
        'biocontainers/nanoq:0.10.0--h031d066_2'}"

    input:
    tuple val(meta), path(ontreads)
    val(output_type) // u: uncompressed; b: Bzip2; g: Gzip; l: Lzma

    output:
    tuple val(meta), path("*.stats")                                                  , emit: stats
    tuple val(meta), path("*")                                                        , emit: reads
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_type = ontreads.baseName.split("\\.")[1]

    if (output_type == 'u') {
        output_ext = ""
    } else if (output_type == 'b') {
        output_ext = ".bz2"
    } else if (output_type == 'g') {
        output_ext = ".gz"
    } else if (output_type == 'l') {
        output_ext = ".lzma"
    }

    """
    nanoq -i $ontreads \\
        -r ${prefix}.stats \\
        -O $output_type \\
        -o ${prefix}_nanoq.${seq_type}${output_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_type = ontreads.baseName.split("\\.")[1]
    if (output_type == 'u') {
        output_ext = ""
    } else if (output_type == 'b') {
        output_ext = ".bz2"
    } else if (output_type == 'g') {
        output_ext = ".gz"
    } else if (output_type == 'l') {
        output_ext = ".lzma"
    }

    """
    touch ${prefix}_nanoq.${seq_type}${output_ext}
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """
}
