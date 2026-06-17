process METAPHLAN3_MERGEMETAPHLANTABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:3.0.12--pyhb7b1952_0' :
        'quay.io/biocontainers/metaphlan:3.0.12--pyhb7b1952_0' }"

    input:
    tuple val(meta), path(profiles)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('metaphlan3'), eval("metaphlan --version 2>&1 | cut -d ' ' -f 3"), emit: versions_metaphlan3, topic: versions


    script:
    def args  = task.ext.args   ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def input = profiles.sort{profile -> profile.toString()}.join(" ")
    """
    merge_metaphlan_tables.py \\
        $args \\
        -o ${prefix}.txt \\
        ${input}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
