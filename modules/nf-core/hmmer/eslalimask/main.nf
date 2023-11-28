process HMMER_ESLALIMASK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1':
        'biocontainers/hmmer:3.3.2--h1b792b2_1' }"

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
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
        $args \\
        $fmask_rfarg \\
        $fmask_allarg \\
        $gmask_rfarg \\
        $gmask_allarg \\
        $pmask_rfarg \\
        $pmask_allarg \\
        -o ${prefix}.masked.sthlm \\
        $unmaskedaln \\
        $maskfile

    gzip ${prefix}.*mask*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer/easel: \$(esl-reformat -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """
}
