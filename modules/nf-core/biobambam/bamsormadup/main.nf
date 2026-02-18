process BIOBAMBAM_BAMSORMADUP {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.185--h85de650_1'
        : 'biocontainers/biobambam:2.0.185--h85de650_1'}"

    input:
    tuple val(meta), path(bams, stageAs: "?/*")
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.bam.bai"), optional: true, emit: bam_index
    tuple val(meta), path("*.cram"), optional: true, emit: cram
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    tuple val("${task.process}"), val('biobambam'), eval("bamsormadup --version |& sed '1!d; s/.*version //; s/.\$//'"), topic: versions, emit: versions_biobambam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("outputformat=cram") ? "cram" : "bam"
    def input_string = bams instanceof List ? bams.join(" I=") : bams
    if (args.contains("outputformat=cram") && fasta == null) {
        error("Reference required for CRAM output.")
    }

    """
    bamcat \\
        I=${input_string} \\
        level=0 \\
    | bamcollate2 \\
        level=0 \\
        ${args2} \\
    | bamsormadup \\
        ${args} \\
        M=${prefix}.metrics.txt \\
        tmpfile=${prefix} \\
        threads=${task.cpus} \\
        > ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("outputformat=cram") ? "cram" : "bam"
    if (args.contains("outputformat=cram") && fasta == null) {
        error("Reference required for CRAM output.")
    }

    """
    touch ${prefix}.${suffix}
    touch ${prefix}.metrics.txt
    """
}
