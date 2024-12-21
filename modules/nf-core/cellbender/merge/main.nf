process CELLBENDER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/055d333b219227f029403ab8591440367f659968eced58b79bdb67c33f30d1e5/data':
        'community.wave.seqera.io/library/cellbender:0.3.0--d6088f41c6eb0a0f' }"

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
