process GCTA_FILTERGRMWITHKEEP {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(grm_id), path(grm_bin), path(grm_n_bin)
    tuple val(meta2), path(keep_file)

    output:
    tuple val(meta), path("*_unrel.grm.id"), path("*_unrel.grm.bin"), path("*_unrel.grm.N.bin"), emit: filtered_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | grep 'version v' | tr -s ' ' | cut -d' ' -f3 | sed 's/^v//'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gcta \\
        --grm ${meta.id} \\
        --keep ${keep_file} \\
        --make-grm \\
        --out ${prefix}_unrel \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unrel.grm.id
    touch ${prefix}_unrel.grm.bin
    touch ${prefix}_unrel.grm.N.bin
    """
}
