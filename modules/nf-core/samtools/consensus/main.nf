process SAMTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

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
