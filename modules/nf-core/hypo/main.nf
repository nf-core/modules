process HYPO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hypo:1.0.3--h9a82719_1':
        'quay.io/biocontainers/hypo:1.0.3--h9a82719_1' }"

    input:
    tuple val(meta), path(sr_bam)
    tuple val(meta2), path(reads)
    path draft
    val genome_size
    val reads_coverage

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val("hypo"), val("1.0.3"), topic: versions, emit: versions_hypo


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${reads.join(' ')}" | tr " " "\\n" > sr.fofn
    hypo \\
        -r @sr.fofn \\
        -d ${draft} \\
        -b ${sr_bam} \\
        -c ${reads_coverage} \\
        -s ${genome_size} \\
        -t ${task.cpus} \\
        -o ${prefix}.fasta \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}
