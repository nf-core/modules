process OMEROBIFROST_PUSH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/omero-py_pip_omero-bifrost:6a0903c81adc5977"

    input:
    tuple val(meta), path(image_folder_path) //stageAs: 'omero_upload/*'
    tuple val(meta2), val(omero_dataset_id)
    tuple val(meta3), path(omero_config_file)

    output:
    tuple val(meta), path("*.xml"), emit: ids_xml
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    omero-bifrost \\
        push \\
        img-folder \\
        ${image_folder_path} \\
        $omero_dataset_id \\
        --to-file \\
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
    touch ./stub_omero_bifrost_output.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omerobifrost: \$(omero-bifrost --version |& sed '1!d ; s/omero-bifrost //')
    END_VERSIONS
    """
}
