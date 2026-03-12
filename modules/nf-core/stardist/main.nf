process STARDIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/stardist:0.9.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.stardist.tif"), emit: mask
    tuple val("${task.process}"), val('stardist'), eval("python -c \"import stardist; print(stardist.__version__)\""), emit: versions_stardist, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('tensorflow'), eval("python -m pip show --version tensorflow | grep \"Version\" | sed -e \"s/Version: //g\""), emit: versions_tensorflow, topic: versions
    tuple val("${task.process}"), val('tifffile'), eval("python -m pip show --version tifffile | grep \"Version\" | sed -e \"s/Version: //g\""), emit: versions_tifffile, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args   ?: ''
    """
    stardist-predict2d \\
        -i $image \\
        -o . \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stardist.tif
    """
}
