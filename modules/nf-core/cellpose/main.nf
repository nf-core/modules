process CELLPOSE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff1ae659d63c642a7ba58bc3c22ca3ee692f5d654ab8eafcfdfbbce010ce1988/data' :
        'community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e' }"

    input:
    tuple val(meta), path(image)
    path(model)

    output:
    tuple val(meta), path("${prefix}/*masks.tif"), emit: mask
    tuple val(meta), path("${prefix}/*flows.tif"), emit: flows, optional: true
    tuple val(meta), path("${prefix}/*seg.npy"), emit: cells, optional: true
    tuple val("${task.process}"), val('cellpose'), eval("pip show cellpose | sed -n 's/^Version: //p'"), topic: versions, emit: versions_cellpose
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('torch'), eval("pip show torch | sed -n 's/^Version: //p'"), topic: versions, emit: versions_torch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model_command = model ? "--pretrained_model ${model}" : ""
    def gpu_flag = task.ext.use_gpu ? "--use_gpu" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    # Container runs as root with HOME=/ which is not writable
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export CELLPOSE_LOCAL_MODELS_PATH=\$PWD/.cellpose
    mkdir -p \$MPLCONFIGDIR \$CELLPOSE_LOCAL_MODELS_PATH

    cellpose \\
        --image_path ${image} \\
        --save_tif \\
        --verbose \\
        ${gpu_flag} \\
        ${model_command} \\
        ${args}

    mkdir -p ${prefix}
    mv *masks.tif ${prefix}/
    mv *flows.tif ${prefix}/ 2>/dev/null || true
    mv *seg.npy ${prefix}/ 2>/dev/null || true
    """

    stub:
    def name = image.name
    def base = name.lastIndexOf('.') != -1 ? name[0..name.lastIndexOf('.') - 1] : name
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/${base}_cp_masks.tif
    """
}
