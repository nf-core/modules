process GCTA_FILTERGRMWITHKEEP {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(grm_id), path(grm_bin), path(grm_n_bin)
    tuple val(meta2), path(keep_file)

    output:
    tuple val(meta), path("${meta.id}_unrel.grm.id"), path("${meta.id}_unrel.grm.bin"), path("${meta.id}_unrel.grm.N.bin"), emit: filtered_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gcta \\
        --grm ${meta.id} \\
        --keep ${keep_file} \\
        --make-grm \\
        --out ${meta.id}_unrel \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_unrel.grm.id
    touch ${meta.id}_unrel.grm.bin
    touch ${meta.id}_unrel.grm.N.bin
    """
}
