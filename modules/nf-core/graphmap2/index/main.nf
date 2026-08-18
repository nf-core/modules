process GRAPHMAP2_INDEX {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'quay.io/biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    path fasta

    output:
    path "*.gmidx", emit: index
    tuple val("${task.process}"), val('graphmap2'), eval("graphmap2 2>&1 | sed -n 's/Version: v//p'"), emit: versions_graphmap2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    graphmap2 \\
        align \\
        -t ${task.cpus} \\
        -I \\
        ${args} \\
        -r ${fasta}
    """

    stub:

    """
    touch ${fasta}.gmidx
    """
}
