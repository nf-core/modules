process CADDSV_GET {
    tag 'CADDSV annotations'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/caddsv:2.0--pyh84cbfca_0':
        'quay.io/biocontainers/caddsv:2.0--pyh84cbfca_0' }"

    output:
    path "caddsv_annotations", emit: annotations
    tuple val("${task.process}"), val('caddsv'), eval("python -c \"import importlib.metadata as m; print(m.version('caddsv'))\""), emit: versions_caddsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    caddsv get \\
        --annotations-dir "caddsv_annotations" \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args

    mkdir -p caddsv_annotations
    touch caddsv_annotations/stub.txt
    """
}
