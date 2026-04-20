process GCTA_REMOVERELATEDSUBJECTS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(grm_id), path(grm_bin), path(grm_n_bin)

    output:
    tuple val(meta), path("*_unrel05.grm.id"), path("*_unrel05.grm.bin"), path("*_unrel05.grm.N.bin"), emit: grm_files
    tuple val(meta), path("*_unrel05.grm.id"), emit: keep_file
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | grep 'version v' | tr -s ' ' | cut -d' ' -f3 | sed 's/^v//'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gcta \\
        --grm ${meta.id} \\
        --grm-cutoff 0.05 \\
        --make-grm \\
        --out ${prefix}_unrel05 \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unrel05.grm.id
    touch ${prefix}_unrel05.grm.bin
    touch ${prefix}_unrel05.grm.N.bin
    """
}
