process AMPLIFY_PREDICT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/amplify:2.0.0--py36hdfd78af_1':
        'biocontainers/amplify:2.0.0--py36hdfd78af_1' }"

    input:
    tuple val(meta), path(faa)
    path(model_dir)

    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    tuple val("${task.process}"), val('AMPlify'), eval( "AMPlify --help 2>&1 | sed -n 's/AMPlify v//p'" ), emit: versions_amplify, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def custom_model_dir = model_dir ? "-md ${model_dir}" : ""
    """
    AMPlify \\
        ${args} \\
        ${custom_model_dir} \\
        -s '${faa}'

    #rename output, because tool includes date and time in name
    mv *.tsv ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
