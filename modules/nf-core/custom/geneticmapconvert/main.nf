process CUSTOM_GENETICMAPCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/062aabd31ebac6f139125e485d5566e928c1b79caf488daa596df02bd1ccbf23/data':
        'community.wave.seqera.io/library/r-data.table_r-janitor_r-r.utils:c8ebef5bb002374e' }"

    input:
    tuple val(meta), path(map_file)

    output:
    tuple val(meta), path("${prefix}.glimpse.map"), emit: glimpse_map
    tuple val(meta), path("${prefix}.plink.map")  , emit: plink_map
    tuple val(meta), path("${prefix}.stitch.map") , emit: stitch_map
    tuple val(meta), path("${prefix}.minimac.map"), emit: minimac_map
    path "versions.yml", emit: versions_geneticmapconvert, topic: versions

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
    END_VERSIONS
    """
}
