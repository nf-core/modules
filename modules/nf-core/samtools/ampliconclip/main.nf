process SAMTOOLS_AMPLICONCLIP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam)
    path bed
    val save_cliprejects
    val save_clipstats

    output:
    tuple val(meta), path("*.clipallowed.bam"), emit: bam
    tuple val(meta), path("*.clipstats.txt"), optional: true, emit: stats
    tuple val(meta), path("*.cliprejects.bam"), optional: true, emit: rejects_bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects = save_cliprejects ? "--rejects-file ${prefix}.cliprejects.bam" : ""
    def stats = save_clipstats ? "-f ${prefix}.clipstats.txt" : ""
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        ampliconclip \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${rejects} \\
        ${stats} \\
        -b ${bed} \\
        -o ${prefix}.clipallowed.bam \\
        ${bam}
    """

    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects = save_cliprejects ? "touch ${prefix}.cliprejects.bam" : ""
    def stats = save_clipstats ? "touch ${prefix}.clipstats.txt" : ""

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    touch ${prefix}.clipallowed.bam
    ${rejects}
    ${stats}
    """
}
