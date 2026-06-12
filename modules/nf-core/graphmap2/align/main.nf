process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'quay.io/biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  index

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val("${task.process}"), val('graphmap2'), eval("graphmap2 2>&1 | sed -n 's/Version: v//p'"), emit: versions_graphmap2, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graphmap2 \\
        align \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        -i ${index} \\
        -d ${reads} \\
        -o ${prefix}.sam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.sam
    """

}
