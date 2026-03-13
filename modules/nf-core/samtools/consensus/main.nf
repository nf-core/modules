process SAMTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752'}"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("*.fastq"), emit: fastq, optional: true
    tuple val(meta), path("*.pileup"), emit: pileup, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-f fastq")
        ? "fastq"
        : args.contains("-f pileup")
            ? "pileup"
            : args.contains("-f fasta")
                ? "fasta"
                : "fasta"

    """
    samtools \\
        consensus \\
        ${args} \\
        -@ ${task.cpus} \\
        -o ${prefix}.${extension} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-f fastq")
        ? "fastq"
        : args.contains("-f pileup")
            ? "pileup"
            : args.contains("-f fasta")
                ? "fasta"
                : "fasta"

    """
    touch ${prefix}.${extension}
    """
}
