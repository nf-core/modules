process ABYSS_ABYSSPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abyss:2.3.10--hf316886_1':
        'biocontainers/abyss:2.3.10--hf316886_1' }"

    input:
    tuple val(meta), path(reads), path(merged)
    val kmersize

    output:
    tuple val(meta), path("*-contigs.fa"),   emit: contigs
    tuple val(meta), path("*-scaffolds.fa"), emit: scaffolds
    tuple val(meta), path("*-stats"),        emit: stats
    tuple val(meta), path("*-abyss.log"),            emit: log
    tuple val("${task.process}"), val('abyss'), eval("abyss-pe version | grep abyss | cut -d\" \" -f3"), topic: versions, emit: versions_abyss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = (task.memory.toGiga()*0.8).toInteger()

    def input_reads = ""
    if (merged.name != 'input.1') {
        input_reads = "in='${reads[0]} ${reads[1]}' se='${merged}'"
    } else {
        input_reads = "in='${reads[0]} ${reads[1]}'"
    }

    """
    abyss-pe \\
        $args \\
        j=$task.cpus \\
        B=${memory}G \\
        k=$kmersize \\
        $input_reads \\
        name=$prefix > ${prefix}-abyss.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}-contigs.fa
    touch ${prefix}-scaffolds.fa
    touch ${prefix}-stats
    touch ${prefix}-abyss.log
    """
}
