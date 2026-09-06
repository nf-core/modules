process STECFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stecfinder:1.1.2--pyhdfd78af_0':
        'quay.io/biocontainers/stecfinder:1.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('stecfinder'), eval("stecfinder --version 2>&1 | sed 's/^.*STECFinder version: //;'"), topic: versions, emit: versions_stecfinder

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stecfinder \\
        -i ${seqs} \\
        ${args} \\
        -t ${task.cpus} > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """

}
