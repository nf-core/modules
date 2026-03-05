process COOLER_CLOAD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contacts), path(index)
    tuple val(meta2), path(chromsizes)
    val(mode)
    val(cool_bin)

    output:
    tuple val(meta), path("*.cool"), emit: cool
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nproc  = mode in ["pairix", "tabix"] ? "--nproc ${task.cpus}" : ""
    """
    cooler cload ${mode} \\
        ${args} \\
        ${nproc} \\
        ${chromsizes}:${cool_bin} \\
        ${contacts} \\
        ${prefix}.cool
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cool
    """
}
