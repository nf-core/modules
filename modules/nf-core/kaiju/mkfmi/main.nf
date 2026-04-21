process KAIJU_MKFMI {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0'
        : 'biocontainers/kaiju:1.10.0--h43eeafb_0'}"

    input:
    tuple val(meta), path(fasta)
    path nodes_dmp, stageAs: "nodes.dmp"
    path names_dmp, stageAs: "names.dmp"
    val keep_intermediate

    output:
    tuple val(meta), path("*.{fmi,dmp}", includeInputs: true), emit: fmi
    tuple val(meta), path("*.bwt"), optional: true, emit: bwt
    tuple val(meta), path("*.sa"), optional: true, emit: sa
    tuple val("${task.process}"), val('kaiju'), eval("kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'"), emit: versions_kaiju, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_cleanup = keep_intermediate ? "" : "rm -f *.{bwt,sa}"
    """
    kaiju-mkbwt \\
        ${args} \\
        -n ${task.cpus} \\
        -o ${prefix} \\
        ${fasta}

    kaiju-mkfmi ${prefix}

    ${run_cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_cleanup = keep_intermediate ? "" : "rm -f *.{bwt,sa}"
    """
    touch ${prefix}.fmi
    touch ${prefix}.bwt
    touch ${prefix}.sa

    ${run_cleanup}
    """
}
