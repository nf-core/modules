process BLOBTK_SNAIL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/08833d1b91f41024e06e2cb5a982598531199c04e6544885d42ef2cb0480de18/data' :
        'community.wave.seqera.io/library/blobtk:0.8.0--2fe0d833a26e0cd9' }"

    input:
    tuple val(meta), path(fasta), path(full_table)
    tuple val(meta2), path(blob_dir)
    val(image_format)

    output:
    tuple val(meta), path("*.{svg,png}"), emit: images
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix              = task.ext.prefix   ?: "${meta.id}"
    def args            = task.ext.args     ?: ""
    def blob_dir_args   = blob_dir          ? "--blobdir ${blob_dir}" : ""
    def full_table_args = full_table        ? "--busco ${full_table}"  : ""
    def output_image    = "${prefix}_snail.${image_format}"

    """
    blobtk snail \\
        --fasta ${fasta} \\
        ${blob_dir_args} \\
        ${full_table_args} \\
        ${args} \\
        --output ${output_image}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svg
    """
}
