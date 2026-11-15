process SVDSS_SEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdss:2.1.1--he17396a_0' :
        'quay.io/biocontainers/svdss:2.1.1--he17396a_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.sfs.txt"), emit: sfs
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, v//'"), emit: versions_svdss, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_flag = reads.name.endsWith('.bam') ? '--bam' : '--fastx'
    """
    SVDSS search \\
        --threads ${task.cpus} \\
        --index ${index} \\
        ${input_flag} ${reads} \\
        ${args} \\
        > ${prefix}.sfs.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.sfs.txt
    """
}
