process OMEROBIFROST_PULL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/omero-py_pip_omero-bifrost:6a0903c81adc5977"

    input:
    tuple val(meta), path(image_id_list)
    tuple val(meta2), path(omero_config_file)

    output:
    tuple val(meta), path("*.ome.tiff"), emit: ometiff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    omero-bifrost \\
        pull \\
        ome-tiffs \\
        --list $image_id_list \\
        --config $omero_config_file \\
        ./ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omerobifrost: \$(omero-bifrost --version |& sed '1!d ; s/omero-bifrost //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ./stub_img.ome.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omerobifrost: \$(omero-bifrost --version |& sed '1!d ; s/omero-bifrost //')
    END_VERSIONS
    """
}
