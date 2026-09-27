process FOLDDISCO_QUERY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/folddisco:2.9375a2d--hb42e459_0':
        'quay.io/biocontainers/folddisco:2.9375a2d--hb42e459_0' }"

    input:
    tuple val(meta), path(pdb), val(query)
    tuple val(meta2), path(index), path(structures, stageAs: "structures/*")

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: results
    tuple val("${task.process}"), val('folddisco'), eval("folddisco version"), topic: versions, emit: versions_folddisco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    folddisco \\
        query \\
        --pdb ${pdb} \\
        --query ${query} \\
        --index ${index}/${index} \\
        --output ${prefix}.tsv \\
        --threads $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
