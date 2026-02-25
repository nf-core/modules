process CELLBENDER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c472f392e7bbfc26dca0088c3bd2349aad10f41b002d1288f6b5958c0951d8df/data':
        'community.wave.seqera.io/library/cellbender_python_webcolors:286b10a91af05a58' }"

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

    """
    echo ${output_layer}
    """

    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        cellbender: \$(cellbender --version)
    END_VERSIONS
    """
}
