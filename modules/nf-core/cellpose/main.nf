process CELLPOSE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "dongzehe/cellpose-gpu:4.0.8"

    input:
    tuple val(meta), path(image)
    path(model)

    output:
    tuple val(meta), path("${prefix}/*masks.tif"), emit: mask
    tuple val(meta), path("${prefix}/*flows.tif"), emit: flows, optional: true
    tuple val(meta), path("${prefix}/*seg.npy"), emit: cells, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CELLPOSE module does not support conda. Please use Docker / Singularity / Podman instead.")
    }
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(cellpose --version | awk 'NR==1 {print \$NF}')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CELLPOSE module does not support conda. Please use Docker / Singularity / Podman instead.")
    }
    def name = image.name
    def base = name.lastIndexOf('.') != -1 ? name[0..name.lastIndexOf('.') - 1] : name
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/${base}_cp_masks.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(cellpose --version | awk 'NR==1 {print \$NF}')
    END_VERSIONS
    """
}
