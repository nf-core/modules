process FQ_GENERATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'quay.io/biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    val meta

    output:
    tuple val(meta), path("*.fastq.gz")                                                                      , emit: fastq
    tuple val("${task.process}"), val('fq'), eval("fq generate --version | sed 's/fq-generate //;s/ .*//'"), emit: versions_fq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq generate \\
        $args \\
        ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_R1.fastq.gz
    echo "" | gzip > ${prefix}_R2.fastq.gz
    """
}
