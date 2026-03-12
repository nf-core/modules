process SAMTOOLS_SORMADUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam")      , emit: bam,  optional: true
    tuple val(meta), path("*.cram")     , emit: cram, optional: true
    tuple val(meta), path("*.csi")      , emit: csi,  optional: true
    tuple val(meta), path("*.crai")     , emit: crai, optional: true
    tuple val(meta), path("*.metrics")  , emit: metrics
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

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
        - \\
    | \\
    samtools markdup \\
        -T ${prefix} \\
        -f ${prefix}.metrics \\
        --threads $task.cpus \\
        $reference \\
        $args5 \\
        - \\
        ${prefix}.${extension}

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
    """
}
