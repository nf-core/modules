
process BIOFORMATS2RAW {
    tag "$meta.id"
    label 'process_high'

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
        --max-workers $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats2raw: \$(bioformats2raw --version |& sed -n '1s/Version = //p')
        bio-formats: \$(bioformats2raw --version |& sed -n '2s/Bio-Formats version = //p')
        ngff: \$(bioformats2raw --version |& sed -n '3s/NGFF specification version = //p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}.ome.zarr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats2raw: \$(bioformats2raw --version |& sed -n '1s/Version = //p')
        bio-formats: \$(bioformats2raw --version |& sed -n '2s/Bio-Formats version = //p')
        ngff: \$(bioformats2raw --version |& sed -n '3s/NGFF specification version = //p')
    END_VERSIONS
    """
}
