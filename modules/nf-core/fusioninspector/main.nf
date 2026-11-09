// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FUSIONINSPECTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bd1e83e196b8989930525832da96ad1535bb36ef04ea5cba53be98c371b58a87/data'
        : 'community.wave.seqera.io/library/fusion-inspector_igv-reports_perl-json-xs_pysam_pruned:9f1b4c0a736e1d8a'}"

    input:
    tuple val(meta), path(reads), path(fusion_list)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*FusionInspector.fusions.tsv"), emit: tsv         , optional:true
    tuple val(meta), path("fi_workdir/*.gtf")            , emit: out_gtf     , optional:true
    tuple val(meta), path("*FusionInspector.log")        , emit: log         , optional:true
    tuple val(meta), path("*html")                       , emit: html        , optional:true
    tuple val(meta), path("*abridged.tsv")               , emit: abridged_tsv, optional:true
    tuple val(meta), path("IGV_inputs")                  , emit: igv_inputs  , optional:true
    tuple val(meta), path("fi_workdir")                  , emit: fi_workdir  , optional:true
    tuple val(meta), path("chckpts_dir")                 , emit: chckpts_dir , optional:true
    tuple val("${task.process}"), val('fusion-inspector'), eval("FusionInspector --version |& sed -n 's/.*version: //p'"), topic: versions, emit: versions_fusioninspector

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta  = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    """
    FusionInspector \\
        --fusions ${fusion_list} \\
        --genome_lib ${reference} \\
        ${fasta} \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix ${prefix} \\
        --vis ${args} ${args2}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.abridged.tsv
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.fusion_inspector_web.html
    mkdir -p chckpts_dir
    mkdir -p fi_workdir
    touch fi_workdir/${prefix}.gtf
    mkdir -p IGV_inputs
    """
}
