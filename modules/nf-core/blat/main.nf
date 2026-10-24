// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BLAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f1c01ef5544da297e213b49ad010a6de4bb152c1bd34881a8b65cb825192ed5/data'
        : 'community.wave.seqera.io/library/ucsc-blat:472--3b87dacc9a173a42' }"

    input:
    tuple val(meta) , path(query)
    tuple val(meta2), path(subject)

    output:
    tuple val(meta), path("*.psl"), emit: psl
    tuple val("${task.process}"), val("blat"), eval("blat 2>&1 | sed '1!d;s/^.*BLAT v. //;s/ fast.*//'"), topic: versions, emit: versions_blat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzip = query.toString().endsWith(".gz")

    """
    in=$query
    if $unzip
    then
        gunzip -cdf $query > ${prefix}.fasta
        in=${prefix}.fasta
    fi

    blat \\
        $args \\
        $subject \\
        \$in \\
        ${prefix}.psl

    if $unzip
    then
        rm ${prefix}.fasta
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.psl
    """
}
