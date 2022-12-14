process ODGI_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::odgi=0.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.8.0--py310hc8f18ef_0':
        'quay.io/biocontainers/odgi:0.8.0--py310hc8f18ef_0' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.og.stats.tsv") , optional: true, emit: tsv
    tuple val(meta), path("*.og.stats.yaml"), optional: true, emit: yaml
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = ".og.stats.tsv"
    if (args.contains("--yaml") || args.contains("--multiqc")) {
        suffix = "og.stats.yaml"
    }
    """
    odgi \\
        stats \\
        --threads $task.cpus \\
        --idx ${graph} \\
        $args > ${prefix}.$suffix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(echo \$(odgi version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
