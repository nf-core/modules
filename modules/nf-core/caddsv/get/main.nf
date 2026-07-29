process CADDSV_GET {
    tag "CADDSV ${flag}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/7580fd97186cafd35f9b898de8ef33a503b42649c99b9b2a50d4cf4eada0bd0d/data'
:         'community.wave.seqera.io/library/caddsv:2.0.1--1bd7ba3bc0ff7a4e' }"

    input:
    val flag

    output:
    path "caddsv_annotations", emit: annotations
    tuple val("${task.process}"), val('caddsv'), eval("caddsv --version"), emit: versions_caddsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (!['annotations', 'segmentnt'].contains(flag)) {
        error("Invalid caddsv get flag: '${flag}'. Expected 'annotations' or 'segmentnt'.")
    }
    """
    caddsv get ${flag} \\
        --annotations-dir "caddsv_annotations" \\
        ${args}
    """

    stub:
    """
    mkdir -p caddsv_annotations
    touch caddsv_annotations/stub.txt
    """
}
