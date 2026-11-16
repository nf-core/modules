process BULLSEYE_GLM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0ab90f0b815b17326d7f16e4ee3399e07bddbff348dc0562f98961140cf2c97/data'
        : 'community.wave.seqera.io/library/r-aod_r-base_r-data.table_r-dplyr_pruned:9ff2262cd35fb547'}"

    input:
    tuple val(meta), path(coverage), path(mutations), path(samplesheet)

    output:
    tuple val(meta), path("*.glm.tsv")            , emit: results
    tuple val(meta), path("*.R_sessionInfo.log")  , emit: session_info
    path "versions.yml"                         , emit: versions_bullseye_glm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'bullseye_glm.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.glm.tsv
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        r-aod: \$(Rscript -e "cat(as.character(packageVersion('aod')))")
    END_VERSIONS
    """
}
