process CADDSV_GET {
    tag "CADDSV ${flag}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/caddsv:2.0--pyh84cbfca_0'
        : 'quay.io/biocontainers/caddsv:2.0--pyh84cbfca_0'}"

    input:
    val flag

    output:
    path "caddsv_annotations", emit: annotations
    tuple val("${task.process}"), val('caddsv'), eval("python -c \"import importlib.metadata as m; print(m.version('caddsv'))\""), emit: versions_caddsv, topic: versions

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
