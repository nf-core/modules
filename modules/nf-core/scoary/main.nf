process SCOARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scoary:1.6.16--py_2' :
        'quay.io/biocontainers/scoary:1.6.16--py_2' }"

    input:
    tuple val(meta), path(genes), path(traits)
    path(tree)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('scoary'), eval('scoary --version 2>&1'), emit: versions_scoary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def newick_tree = tree ? "-n ${tree}" : ""
    """
    scoary \\
        ${args} \\
        --no-time \\
        --threads ${task.cpus} \\
        --traits ${traits} \\
        --genes ${genes} \\
        ${newick_tree}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
