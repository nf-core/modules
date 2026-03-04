process STARDIST {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d964e0bef867bb2ff1a309c9c087d8d83ac734ce3aa315dd8311d4c1bfdafd8e/data' :
        'community.wave.seqera.io/library/python_pip_imagecodecs_nvidia-cublas-cu12_pruned:b668bcb6d531d350' }"

    input:
    tuple val(meta), path(image)
    tuple val(model_name), path(model_path)

    output:
    tuple val(meta), path("*.stardist.tif"), emit: mask
    tuple val("${task.process}"), val('stardist'), eval("pip show stardist | grep '^Version:' | sed 's/Version: //'"), topic: versions, emit: versions_stardist
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tensorflow'), eval("pip show tensorflow | grep '^Version:' | sed 's/Version: //'"), topic: versions, emit: versions_tensorflow
    tuple val("${task.process}"), val('tifffile'), eval("pip show tifffile | grep '^Version:' | sed 's/Version: //'"), topic: versions, emit: versions_tifffile

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model_command = model_path ? "-m ${model_path}" : model_name ? "-m ${model_name}" : ""
    """
    stardist-predict2d \\
        -i $image \\
        -o . \\
        $model_command \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stardist.tif
    """
}
