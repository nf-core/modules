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
        r-nfcore.utils: \$(Rscript -e "cat(as.character(packageVersion('nfcore.utils')))")
        r-gprofiler2: \$(Rscript -e "cat(as.character(packageVersion('gprofiler2')))")
        gprofiler-data: \$(Rscript -e "cat(gprofiler2::get_version_info()[['gprofiler_version']])")
    END_VERSIONS
    """
}
