process CENTRIFUGER_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuger:1.1.0--hf426362_0':
        'biocontainers/centrifuger:1.1.0--hf426362_0' }"

    input:    
    tuple val(meta), path(reference)
    path conversion_table
    path taxonomy_nodes
    path taxonomy_names


    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.*"), emit: db
    tuple val("${task.process}"), val("centrifuger"), eval("centrifuger -v 2>&1 | head -n 1"), emit: versions_centrifuger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

   """

    centrifuger-build \\
        $args \\
        -t ${task.cpus} \\
        -r $reference \\
        --conversion-table $conversion_table \\
        --taxonomy-tree $taxonomy_nodes \\
        --name-table $taxonomy_names \\
        -o ${prefix} 



    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.1.cfr
    touch ${prefix}.2.cfr
    touch ${prefix}.3.cfr
    touch ${prefix}.4.cfr

    """


}
