process GPROFILER2_GCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbb3149e57b187cd63966a1e25bc71a07f9b965509f5c1ffcb81e071bc8a821a/data':
        'community.wave.seqera.io/library/r-gprofiler2_r-nfcore.utils:b62ba45166a0c1ea' }"

    input:
    tuple val(meta), path(ids_tsv), val(target)

    output:
    tuple val(meta), path("*.gprofiler2.gconvert.tsv"), emit: converted_ids
    tuple val(meta), path("*R_sessionInfo.log")        , emit: session_info
    tuple val("${task.process}"), val("r-base")        , eval("Rscript -e 'cat(as.character(getRversion()))'")                                     , topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val("r-nfcore.utils"), eval("Rscript -e 'cat(as.character(packageVersion(\"nfcore.utils\")))'")                   , topic: versions, emit: versions_nfcoreutils
    tuple val("${task.process}"), val("r-gprofiler2")  , eval("Rscript -e 'cat(as.character(packageVersion(\"gprofiler2\")))'")                     , topic: versions, emit: versions_gprofiler2
    tuple val("${task.process}"), val("gprofiler-data"), eval("Rscript -e 'cat(gprofiler2::get_version_info()[[\"gprofiler_version\"]])'")          , topic: versions, emit: versions_gprofilerdata

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'gprofiler2_gconvert.R'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.gprofiler2.gconvert.tsv

    Rscript -e "nfcore.utils::create_log_session_info(
        log_path = '${prefix}.R_sessionInfo.log'
    )"
    """
}
