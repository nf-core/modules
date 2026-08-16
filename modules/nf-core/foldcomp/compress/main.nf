process FOLDCOMP_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldcomp:1.0.0--h7f5d12c_0':
        'quay.io/biocontainers/foldcomp:1.0.0--h7f5d12c_0' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("*fcz")                                                            , emit: fcz
    tuple val("${task.process}"), val('foldcomp'), eval("foldcomp --version | cut -d' ' -f2"), emit: versions_foldcomp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = pdb.isDirectory() ? "_fcz" : ".fcz"
    """
    foldcomp \\
        compress \\
        -t ${task.cpus} \\
        ${pdb} \\
        ${prefix}${ext} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fcz
    """
}
