process SAMTOOLS_SORMADUP {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam")      , emit: bam
    tuple val(meta), path("*.metrics")  , emit: metrics
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def sort_memory = task.memory.mega ?: "4G"

    """
    samtools cat \\
        --threads $task.cpus \\
        ${input}  \\
    | \\
    samtools collate \\
        -O \\
        -u \\
        --threads $task.cpus \\
        ${reference} \\
        - \\
    | \\
    samtools fixmate \\
        -m \\
        -u \\
        --threads $task.cpus \\
        - \\
        - \\
    | \\
    samtools sort \\
        -u \\
        -T ${prefix} \\
        --threads $task.cpus \\
        - \\
    | \\
    samtools markdup \\
        -T ${prefix} \\
        -f ${prefix}.metrics \\
        --threads $task.cpus \\
        $args \\
        - \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
