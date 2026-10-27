// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HMMER_ESLSFETCHINDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41f71771ca71c8a250d63b4557eb1ce70efc3fc28e565297a9f1890e687561cc/data'
:         'community.wave.seqera.io/library/hmmer_coreutils:e129de528776aa4c' }"

    input:
    tuple val(meta), path(seqfile)

    output:
    tuple val(meta), path("${seqfile}.ssi"), emit: ssi
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('easel'), eval("esl-sfetch -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions
    tuple val("${task.process}"), val('coreutils'), eval("sort --version |& sed '1!d ; s/sort (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    esl-sfetch --index $seqfile
    """

    stub:
    """
    touch ${seqfile}.ssi
    """
}
