process CELLPHENO_NIS {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    // CellPheno NIS is a custom LibTorch + OpenCV + CUDA executable built from the
    // source repository (https://github.com/Chrisa142857/Lightsheet_microscopy_image_3D_nuclei_instance_segmentation).
    // It has no Conda/Bioconda package and requires a CUDA runtime, so it is
    // distributed only as a dedicated GPU container image (cf. the parabricks and
    // numorph/3dunet modules, which likewise ship a vendor GPU image with no conda).
    container "quay.io/wzq10101/cellpheno-nis:1.0.0"

    input:
    tuple val(meta), path(tile_dir)
    path models

    output:
    tuple val(meta), path("*_NIScpp_results_*.zip"), emit: nis
    tuple val(meta), path("*_remap.zip")           , emit: remap, optional: true
    tuple val("${task.process}"), val('cellpheno_nis'), eval("cat /usr/local/share/cellpheno-nis/VERSION"), topic: versions, emit: versions_cellpheno_nis

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba: NIS is a
    // GPU binary with no conda package (same guard as the nf-core parabricks modules).
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("The CELLPHENO_NIS module does not support Conda/Mamba. Please use Docker / Singularity / Apptainer / Podman instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The GPU is selected by the executor (e.g. CUDA_VISIBLE_DEVICES); NIS defaults
    // to `cuda:0`. Override with `--device cuda:N` via `task.ext.args` if needed.
    """
    main \\
        --model_root ${models} \\
        --data_root ${tile_dir} \\
        --save_root . \\
        --brain_tag ${prefix} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_NIScpp_results_zmin0_seg_meta.zip
    touch ${prefix}_NIScpp_results_zmin0_binary_mask.zip
    touch ${prefix}_NIScpp_results_zmin0_instance_center.zip
    touch ${prefix}_NIScpp_results_zmin0_instance_coordinate.zip
    touch ${prefix}_NIScpp_results_zmin0_instance_label.zip
    touch ${prefix}_NIScpp_results_zmin0_instance_volume.zip
    touch ${prefix}_remap.zip
    """
}
