process GCTA_ADDGRMS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(mgrm_file), path(grm_files)

    output:
    tuple val(meta), path("${meta.id}.grm.id"), path("${meta.id}.grm.bin"), path("${meta.id}.grm.N.bin"), emit: combined_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gcta \\
        --mgrm ${mgrm_file} \\
        --make-grm \\
        --out ${meta.id} \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    """
    touch ${meta.id}.grm.id
    touch ${meta.id}.grm.bin
    touch ${meta.id}.grm.N.bin
    """
}
