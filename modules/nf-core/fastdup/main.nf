process FASTDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c55070589353b3e1837ca3414c4f182d3674cbf55a64edee07e8bf75370762a9/data':
        'community.wave.seqera.io/library/fastdup:1.0.0--a9b28abff06bb2bb' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('fastdup'), eval("fastdup --version"), topic: versions, emit: versions_fastdup


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${reads}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    fastdup \\
        $args \\
        --input $reads \\
        --metrics ${prefix}.metrics.txt \\
        --output ${prefix}.bam \\
        --num-threads $task.cpus

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_command = args.contains("--index-format CSI") ? "touch ${prefix}.csi"
                        : args.contains("--create-index")   ? "touch ${prefix}.bai" : ""

    """

    touch ${prefix}.bam
    ${index_command}
    touch ${prefix}.metrics.txt

    """
}
