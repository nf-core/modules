process MODKIT_BEDMETHYLTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0':
        'biocontainers/ont-modkit:0.6.1--hcdda2d0_0' }"

    input:
    tuple val(meta),  path(bedmethyl)
    tuple val(meta2), path(chromsizes)
    val modcodes

    output:
    tuple val(meta), path("*.bw"), emit: bw
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def mods      = modcodes instanceof List ? modcodes.join(',') : modcodes
    def input_cmd = bedmethyl.getName().endsWith('.gz') ? "gzip -cd" : "cat"
    """
    $input_cmd $bedmethyl |\\
        modkit bedmethyl tobigwig \\
        $args \\
        --nthreads $task.cpus \\
        --sizes $chromsizes \\
        --mod-codes $mods \\
        - \\
        ${prefix}.bw
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    """
}
