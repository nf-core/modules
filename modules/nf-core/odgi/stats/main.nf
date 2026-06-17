process ODGI_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.9.0--py312h5e9d817_1':
        'quay.io/biocontainers/odgi:0.9.0--py312h5e9d817_1' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.og.stats.tsv") , optional: true, emit: tsv
    tuple val(meta), path("*.og.stats.yaml"), optional: true, emit: yaml
    tuple val("${task.process}"), val('odgi'), eval("odgi version | sed 's/^v//; s/-.*//'"), emit: versions_odgi, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "og.stats.tsv"
    if (args.contains("--yaml") || args.contains("--multiqc")) {
        suffix = "og.stats.yaml"
    }
    """
    odgi \\
        stats \\
        --threads $task.cpus \\
        --idx ${graph} \\
        $args > ${prefix}.$suffix
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.og.stats.tsv
    """
}
