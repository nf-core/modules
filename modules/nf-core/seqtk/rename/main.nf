process SEQTK_RENAME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.gz")     , emit: sequences
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fasta"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        extension = "fastq"
    }
    """
    seqtk \\
        rename \\
        $args \\
        $sequences \\
        $prefix | \\
        gzip -c --no-name > ${prefix}.renamed.${extension}.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = "fasta"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        extension = "fastq"
    }
    """
    echo "" | gzip > ${prefix}.renamed.${extension}.gz
    """
}
