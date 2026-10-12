process HMMER_ESLSFETCHINDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2c9b2c2ded44fc2076629da06a841d3bd8a9f847f6d6e4b6a3aac38897870bc/data' :
        'community.wave.seqera.io/library/coreutils_hmmer:4d607118d1d3ddb5' }"

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
