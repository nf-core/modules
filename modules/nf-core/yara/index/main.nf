process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yara:1.0.2--2' :
        'quay.io/biocontainers/yara:1.0.2--2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}*")   , emit: index
    tuple val("${task.process}"), val('yara'), eval("yara_indexer --version 2>&1 | grep 'yara_indexer version' | sed 's/^.*yara_indexer version: //; s/ .*\$//'"), topic: versions, emit: versions_yara


    script:
    def args   = task.ext.args   ?: ''

    """
    yara_indexer \\
        $fasta \\
        -o ${fasta} \\
        ${args}
    """

    stub:
    """
    touch ${fasta}.lf.{drp,drs,drv,pst}
    touch ${fasta}.rid.{concat,limits}
    touch ${fasta}.sa.{ind,len,val}
    touch ${fasta}.txt.{concat,limits,size}
    """
}
