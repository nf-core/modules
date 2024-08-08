process STARDIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/stardist:0.9.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.stardist.tif"), emit: mask
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args   ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def VERSION            = "0.9.1"
    def VERSION_PYTHON     = "3.9"
    def VERSION_TENSORFLOW = "2.10.0"
    def VERSION_TIFFFILE   = "<2022.4.22"


    """
    stardist-predict2d \\
        -i $image \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: $VERSION
        python: $VERSION_PYTHON
        tensorflow: $VERSION_TENSORFLOW
        tifffile: $VERSION_TIFFFILE
    END_VERSIONS
    """

    stub:
    def args               = task.ext.args   ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def VERSION            = "0.9.1"
    def VERSION_PYTHON     = "3.9"
    def VERSION_TENSORFLOW = "2.10.0"
    def VERSION_TIFFFILE   = "<2022.4.22"
    """
    touch ${prefix}.stardist.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: $VERSION
        python: $VERSION_PYTHON
        tensorflow: $VERSION_TENSORFLOW
        tifffile: $VERSION_TIFFFILE
    END_VERSIONS
    """
}
