// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc124a4b28a7b8cf45e75edd574368c34ac97d966a9f493120b51fe4338536b8/data'
:         'community.wave.seqera.io/library/yara:1.0.2--6d4e380b8f30498b' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}*")   , emit: index
    tuple val("${task.process}"), val('yara'), eval("yara_indexer --version 2>&1 | grep 'yara_indexer version' | sed 's/^.*yara_indexer version: //; s/ .*\$//'"), topic: versions, emit: versions_yara

    when:
    task.ext.when == null || task.ext.when

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
