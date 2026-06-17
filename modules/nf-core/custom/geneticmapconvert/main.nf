process CUSTOM_GENETICMAPCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b3/b314d622eb495f740bcfa163ceda8f65cb9c77afa069f533d41c810c33eeecdb/data':
        'community.wave.seqera.io/library/r-data.table_r-janitor_r-nfcore.utils_r-r.utils:27d6c57bc9981f04' }"

    input:
    tuple val(meta), path(map_file)

    output:
    tuple val(meta), path("${prefix}.glimpse.map")      , emit: glimpse_map
    tuple val(meta), path("${prefix}.plink.map")        , emit: plink_map
    tuple val(meta), path("${prefix}.stitch.map")       , emit: stitch_map
    tuple val(meta), path("${prefix}.minimac.map")      , emit: minimac_map
    tuple val(meta), path("${prefix}.R_sessionInfo.log"), emit: session_info
    path "versions.yml", emit: versions, topic: versions


    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ''

    """
    echo ${args}
    """

    template 'geneticmapconvert.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.glimpse.map
    touch ${prefix}.plink.map
    touch ${prefix}.stitch.map
    touch ${prefix}.minimac.map
    touch ${prefix}.eagle.map

    Rscript -e "nfcore.utils::process_end(
        packages = list('r-data.table' = 'data.table', 'r-janitor' = 'janitor'),
        task_name = '${task.process}',
        versions_path = 'versions.yml',
        log_path = '${prefix}.R_sessionInfo.log'
    )"
    """
}
