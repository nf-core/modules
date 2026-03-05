process SEQFU_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqfu:1.22.3--hfd12232_2'
        : 'biocontainers/seqfu:1.22.3--hfd12232_2'}"

    input:
    // stats can get one or more fasta or fastq files
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    tuple val("${task.process}"), val('seqfu'), eval('seqfu version'), emit: versions_seqfu, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqfu \\
        stats \\
        ${args} \\
        --multiqc ${prefix}_mqc.txt \\
        ${files} > ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.tsv
    touch ${prefix}_mqc.txt
    """
}
