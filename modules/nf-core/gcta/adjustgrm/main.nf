process GCTA_ADJUSTGRM {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(grm_id), path(grm_bin), path(grm_n_bin)

    output:
    tuple val(meta), path("${meta.id}_adj.grm.id"), path("${meta.id}_adj.grm.bin"), path("${meta.id}_adj.grm.N.bin"), emit: grm_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | grep 'version v' | tr -s ' ' | cut -d' ' -f3 | sed 's/^v//'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gcta \\
        --grm ${meta.id} \\
        --grm-adj 0 \\
        --make-grm \\
        --out ${meta.id}_adj \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}_adj.grm.id
    touch ${meta.id}_adj.grm.bin
    touch ${meta.id}_adj.grm.N.bin
    """
}
