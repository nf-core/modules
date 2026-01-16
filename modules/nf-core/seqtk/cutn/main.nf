process SEQTK_CUTN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bed")    , emit: bed
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqtk \\
        cutN \\
        $args \\
        -g $fasta \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """

}
