process SPLITUBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/splitubam:0.1.1--hc9368f3_0':
        'biocontainers/splitubam:0.1.1--hc9368f3_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('splitubam'), eval("splitubam --version | sed 's/.* //'"), emit: versions_splitubam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    splitubam \\
        $args \\
        --threads $task.cpus \\
        $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def match = (args =~ /--split\s+(\d+)/)
    def create_cmd = ""
    if (match) {
        def n_splits = match[0][1].toInteger()
        create_cmd = (1..n_splits).collect { i ->
            def formattedIteration = String.format('%03d', i)
            "touch ${formattedIteration}.${bam}.bam"
        }.join(" ")
    } else { error("No `--split N` detected in args") }
    """
    $create_cmd
    """
}
