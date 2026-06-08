process CUSTOM_GENETICMAPCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2df9ef33f7170686ac09d8c182f46b3bb10fb5d881a7253c530e8bac09337b50/data':
        'community.wave.seqera.io/library/r-data.table_r-janitor_r-nfcore.utils_r-r.utils:458c7966e529a0f9' }"

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed '1!d; s/.*version //; s/ .*//')
        r-data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))")
        r-janitor: \$(Rscript -e "cat(as.character(packageVersion('janitor')))")
        r-nfcore.utils: \$(Rscript -e "cat(as.character(packageVersion('nfcore.utils')))")
    END_VERSIONS
    """
}
