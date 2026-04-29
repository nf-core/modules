process ABYSS_ABYSSPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abyss:2.3.10--hf316886_1':
        'quay.io/biocontainers/abyss:2.3.10--hf316886_1' }"

    input:
    tuple val(meta), path(reads), path(merged)
    val kmersize

    output:
    tuple val(meta), path("*-contigs.fa.gz"),   emit: contigs
    tuple val(meta), path("*-scaffolds.fa.gz"), emit: scaffolds
    tuple val(meta), path("*-stats"),           emit: stats
    tuple val(meta), path("*-abyss.log"),       emit: log
    tuple val("${task.process}"), val('abyss'), eval("abyss-pe version | grep abyss | cut -d\" \" -f3"), topic: versions, emit: versions_abyss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = (task.memory.toGiga()*0.8).toInteger()

    def input_reads = ""
    if (merged) {
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

    if [ -f ${prefix}-contigs.fa ]; then
        gzip -cn ${prefix}-contigs.fa > ${prefix}-contigs.fa.gz
    fi
    if [ -f ${prefix}-scaffolds.fa ]; then
        gzip -cn ${prefix}-scaffolds.fa > ${prefix}-scaffolds.fa.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip -n > ${prefix}-contigs.fa.gz
    echo "" | gzip -n > ${prefix}-scaffolds.fa.gz
    touch ${prefix}-stats
    touch ${prefix}-abyss.log
    """
}
