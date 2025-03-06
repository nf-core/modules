process SAMTOOLS_SORMADUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam")      , emit: bam,  optional: true
    tuple val(meta), path("*.cram")     , emit: cram, optional: true
    tuple val(meta), path("*.csi")      , emit: csi,  optional: true
    tuple val(meta), path("*.crai")     , emit: crai, optional: true
    tuple val(meta), path("*.metrics")  , emit: metrics
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args5.contains("--output-fmt sam") ? "sam" :
                    args5.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    def sort_memory = (task.memory.mega/task.cpus*0.75).intValue()

    """
    samtools cat \\
        $args \\
        ${input}  \\
    | \\
    samtools collate \\
        $args2 \\
        -O \\
        -u \\
        -T ${prefix}.collate \\
        --threads $task.cpus \\
        ${reference} \\
        - \\
    | \\
    samtools fixmate \\
        $args3 \\
        -m \\
        -u \\
        --threads $task.cpus \\
        - \\
        - \\
    | \\
    samtools sort \\
        $args4 \\
        -u \\
        -T ${prefix}.sort \\
        --threads $task.cpus \\
        -m ${sort_memory}M \\
        - \\
    | \\
    samtools markdup \\
        -T ${prefix} \\
        -f ${prefix}.metrics \\
        --threads $task.cpus \\
        $args5 \\
        - \\
        ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args5 = task.ext.args5 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args5.contains("--output-fmt sam") ? "sam" :
                    args5.contains("--output-fmt cram") ? "cram" :
                    "bam"

    """
    touch ${prefix}.${extension}
    touch ${prefix}.metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
