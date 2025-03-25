process IMAGESEGMENTATIONMETRICS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/fire_monai_rasterio_tifffile_pruned:106d8ac1dc78149c'

    input:
    tuple val(meta), path(ref_img), path(images2compare), val(metricname)

    output:
    tuple val(meta), path("*.csv"), emit: segmentation_metrics
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    segmentation_metrics.py run --path_ref ${ref_img} --path_list_seg "${images2compare}" --list_metrics "${metricname}" ${args}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seg-metrics: \$(segmentation_metrics.py version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seg-metrics: \$(segmentation_metrics.py version)

    END_VERSIONS
    """
}
