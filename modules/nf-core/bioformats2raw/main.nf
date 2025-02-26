
process BIOFORMATS2RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/bioformats2raw:0.9.4--3eec45888b3759e5"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*ome.zarr"),         emit: omezarr
    path "versions.yml"               ,         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bioformats2raw \\
        $image \\
        ${prefix}.ome.zarr \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats2raw: \$(bioformats2raw --version |& sed '1!d ; s/bioformats2raw //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.ome.zarr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats2raw: \$(bioformats2raw --version |& sed '1!d ; s/bioformats2raw //')
    END_VERSIONS
    """
}
