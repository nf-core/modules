process MERYL_UNIONSUM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_0':
        'biocontainers/meryl:1.4.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(meryl_dbs)
    val kvalue

    output:
    tuple val(meta), path("*.unionsum.meryl"), emit: meryl_db
    tuple val("${task.process}"), val('meryl'), eval("meryl --version |& sed 's/meryl //'"), emit: versions_meryl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meryl union-sum \\
        k=$kvalue \\
        threads=$task.cpus \\
        memory=${task.memory.toGiga()} \\
        $args \\
        output ${prefix}.unionsum.meryl \\
        $meryl_dbs
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unionsum.meryl
    """
}
