process CELLBENDER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/174a866a84c6a99dddb058d9b8594db83d5dc5e24332b68620f0a1ffa1bfe82f/data':
        'community.wave.seqera.io/library/cellbender_pyyaml:9b97c45868d30e42' }"

    input:
    tuple val(meta), path(filtered), path(unfiltered), path(cellbender_h5)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    output_layer = task.ext.output_layer ?: "cellbender"
    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(cellbender --version)
    END_VERSIONS
    """
}
