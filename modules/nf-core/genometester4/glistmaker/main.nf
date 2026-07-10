process GENOMETESTER4_GLISTMAKER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/36ce7645eda890c172cecfa37b129d3600e0a23f64bbf05465df4423212ac958/data':
        'community.wave.seqera.io/library/genometester4:4.0--061bc5822e0dd41d' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.list"), emit: list
    tuple val("${task.process}"), val('genometester4'), eval("glistmaker --version | sed -e 's/glistmaker v//g' "), topic: versions, emit: versions_genometester4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    glistmaker \\
        $fasta \\
        $args \\
        --num_threads $task.cpus \\
        -o ${prefix}
    mv ${prefix}_22.list ${prefix}.list
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.list
    """
}
