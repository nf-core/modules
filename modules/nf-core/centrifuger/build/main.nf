process CENTRIFUGER_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:    
    tuple val(meta), path(reference)
    path conversion_table
    path taxonomy_nodes
    path taxonomy_names


    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.*"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

   """
    centrifuger-build \\
        -t ${task.cpus} \\
        -r $reference \\
        --conversion-table $conversion_table \\
        --taxonomy-tree $taxonomy_nodes \\
        --name-table $taxonomy_names \\
        -o ${prefix} 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuger: \$(centrifuger-build --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.cfr
    touch ${prefix}.2.cfr
    touch ${prefix}.3.cfr
    touch ${prefix}.4.cfr 
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuger: "stub"
    END_VERSIONS
    """


}
