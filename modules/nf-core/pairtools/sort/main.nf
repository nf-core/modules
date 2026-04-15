process PAIRTOOLS_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.3--py39h7a39fba_0' :
        'biocontainers/pairtools:1.1.3--py39h7a39fba_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: sorted
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*pairtools.*version //'") , emit: versions_pairtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer = task.memory.toGiga().intdiv(2)
    """
    pairtools \\
        sort \\
        $args \\
        --nproc $task.cpus \\
        --memory ${buffer}G \\
        -o ${prefix}.pairs.gz \\
        $input
    """
}
