process VG_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.73.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.73.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.xg")       , emit: xg
    tuple val(meta), path("*.vgi")      , emit: vg_index, optional: true
    tuple val("${task.process}"), val('vg'), eval("vg 2>&1 | sed -n 's/.*version v\\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_vg


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vg index \\
        --temp-dir . \\
        --threads ${task.cpus} \\
        --xg-name ${prefix}.xg \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def vg_index = args.contains('--index-sorted-vg') ? "touch ${prefix}.vg.vgi" : ""

    """
    touch ${prefix}.xg
    ${vg_index}
    """
}
