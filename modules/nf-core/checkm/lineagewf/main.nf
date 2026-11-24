// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CHECKM_LINEAGEWF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75911411b834b3cbd473babea3721a6f8ca3d5c55783cc02cffadf583e8d6dd5/data'
        : 'community.wave.seqera.io/library/checkm-genome:1.4.0--c2e5d2423bcf948a'}"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")
    val fasta_ext
    path db

    output:
    tuple val(meta), path("${prefix}"), emit: checkm_output
    tuple val(meta), path("${prefix}/lineage.ms.gz"), emit: marker_file
    tuple val(meta), path("${prefix}.tsv"), emit: checkm_tsv
    tuple val("${task.process}"), val('checkm'), eval("checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//'"), emit: versions_checkm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def checkm_db = db ? "export CHECKM_DATA_PATH=${db}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${checkm_db}

    checkm \\
        lineage_wf \\
        -t ${task.cpus} \\
        -x ${fasta_ext} \\
        ${args} \\
        input_bins/ \\
        ${prefix}

    cp ${prefix}/quality_report.tsv ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p checkm_dummy_db
    export CHECKM_DATA_PATH=\$PWD/checkm_dummy_db
    mkdir ${prefix}/
    echo "" | gzip > ${prefix}/lineage.ms.gz
    touch ${prefix}/quality_report.tsv ${prefix}.tsv
    """
}
