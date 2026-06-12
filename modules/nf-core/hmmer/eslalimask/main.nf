process HMMER_ESLALIMASK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(unmaskedaln), val(fmask_rf), val(fmask_all), val(gmask_rf), val(gmask_all), val(pmask_rf), val(pmask_all)
    path  maskfile

    output:
    tuple val(meta), path("*.masked.sthlm.gz"), emit: maskedaln
    path  "*.fmask-rf.gz"                     , emit: fmask_rf , optional: true
    path  "*.fmask-all.gz"                    , emit: fmask_all, optional: true
    path  "*.gmask-rf.gz"                     , emit: gmask_rf , optional: true
    path  "*.gmask-all.gz"                    , emit: gmask_all, optional: true
    path  "*.pmask-rf.gz"                     , emit: pmask_rf , optional: true
    path  "*.pmask-all.gz"                    , emit: pmask_all, optional: true
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions
    tuple val("${task.process}"), val('easel'), eval("esl-alimask -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions


    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fmask_rfarg  = fmask_rf  ? "--fmask-rf ${prefix}.fmask-rf"   : ""
    def fmask_allarg = fmask_all ? "--fmask-all ${prefix}.fmask-all" : ""
    def gmask_rfarg  = gmask_rf  ? "--gmask-rf ${prefix}.gmask-rf"   : ""
    def gmask_allarg = gmask_all ? "--gmask-all ${prefix}.gmask-all" : ""
    def pmask_rfarg  = pmask_rf  ? "--pmask-rf ${prefix}.pmask-rf"   : ""
    def pmask_allarg = pmask_all ? "--pmask-all ${prefix}.pmask-all" : ""
    """
    esl-alimask \\
        $fmask_rfarg \\
        $fmask_allarg \\
        $gmask_rfarg \\
        $gmask_allarg \\
        $pmask_rfarg \\
        $pmask_allarg \\
        -o ${prefix}.masked.sthlm \\
        $args $unmaskedaln $maskfile

    gzip ${prefix}.*mask*
    """

    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def fmask_rfarg  = fmask_rf  ? "touch ${prefix}.fmask-rf"   : ""
    def fmask_allarg = fmask_all ? "touch ${prefix}.fmask-all" : ""
    def gmask_rfarg  = gmask_rf  ? "touch ${prefix}.gmask-rf"   : ""
    def gmask_allarg = gmask_all ? "touch ${prefix}.gmask-all" : ""
    def pmask_rfarg  = pmask_rf  ? "touch ${prefix}.pmask-rf"   : ""
    def pmask_allarg = pmask_all ? "touch ${prefix}.pmask-all" : ""

    """
    touch ${prefix}.masked.sthlm
    ${fmask_rfarg}
    ${fmask_allarg}
    ${gmask_rfarg}
    ${gmask_allarg}
    ${pmask_rfarg}
    ${pmask_allarg}

    gzip ${prefix}.*mask*
    """
}
