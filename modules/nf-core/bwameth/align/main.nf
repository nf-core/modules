process BWAMETH_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.9--pyh7e72e81_0' :
        'biocontainers/bwameth:0.2.9--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('bwameth'), eval("bwameth.py --version | cut -f2 -d ' '"), emit: versions_bwameth, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def args2      = task.ext.args2 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    export BWA_METH_SKIP_TIME_CHECKS=1
    ln -sf \$(readlink ${fasta}) ${index}/${fasta}

    bwameth.py \\
        ${args} \\
        ${read_group} \\
        -t ${task.cpus} \\
        --reference ${index}/${fasta} \\
        ${reads} \\
        | samtools view ${args2} -@ ${task.cpus} -bhS -o ${prefix}.bam -
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
