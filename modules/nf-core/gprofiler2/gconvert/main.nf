process GPROFILER2_GCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e4b0e10a72db4ad519c128c4e3cef6e10bc1a83440af31f105ab389a5532589a/data':
        'community.wave.seqera.io/library/r-ggplot2_r-gprofiler2:fab855ea9f680400' }"

    input:
    tuple val(meta), path(ids_tsv), val(target)

    output:
    tuple val(meta), path("*.gprofiler2.gconvert.tsv"), emit: converted_idsg
    tuple val(meta), path("*R_sessionInfo.log")        , emit: session_info
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'gprofiler2_gconvert.R'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.gprofiler2.gconvert.tsv
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-gprofiler2: \$(Rscript -e "cat(as.character(packageVersion('gprofiler2')))")
        gprofiler-data: \$(Rscript -e "cat(gprofiler2::get_version_info()[['gprofiler_version']])")
    END_VERSIONS
    """
}
