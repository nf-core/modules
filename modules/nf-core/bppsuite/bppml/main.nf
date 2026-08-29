process BPPSUITE_BPPML {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bppsuite:3.0.0--hd63eeec_0':
        'quay.io/biocontainers/bppsuite:3.0.0--hd63eeec_0' }"

    input:
    tuple val(meta), path(alignment), path(tree), path(parameters)

    output:
    tuple val(meta), path("*.params*"), emit: parameters
    tuple val(meta), path("*.{nwk,tree,dnd}*"), emit: tree, optional: true
    tuple val(meta), path("*.backup*"), emit: backup, optional: true
    tuple val(meta), path("*.infos*"), emit: infos, optional: true
    tuple val(meta), path("*.messages*"), emit: messages, optional: true
    tuple val(meta), path("*.profile*"), emit: profile, optional: true
    tuple val("${task.process}"), val('bppml'), eval("bppml | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_bppsuite

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bppml \\
        NUMBER_OF_THREADS=${task.cpus} \\
        ALIGNMENT=${alignment} \\
        TREE=${tree} \\
        DATA=${prefix} \\
        param=${parameters} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.params.txt
    touch ${prefix}.ML.dnd_1
    touch ${prefix}.infos_1
    touch ${prefix}.messages
    touch ${prefix}.profile
    """
}
