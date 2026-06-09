process CUSTOM_GENETICMAPCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/af/afe3fd20fa8b3096a9fc6b9b2725d92516a036651545a5f9ec1e6a53b5a365d9/data':
        'community.wave.seqera.io/library/r-data.table_r-janitor_r-nfcore.utils_r-r.utils:e881939d5868e7db' }"

    input:
    tuple val(meta), path(map_file)

    output:
    tuple val(meta), path("${prefix}.glimpse.map"), emit: glimpse_map
    tuple val(meta), path("${prefix}.plink.map")  , emit: plink_map
    tuple val(meta), path("${prefix}.stitch.map") , emit: stitch_map
    tuple val(meta), path("${prefix}.minimac.map"), emit: minimac_map
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

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

    Rscript -e "nfcore.utils::process_end(packages = list('r-data.table' = 'data.table', 'r-janitor' = 'janitor'), task_name = '${task.process}')"
    """
}
