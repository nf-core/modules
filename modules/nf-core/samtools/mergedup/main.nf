process SAMTOOLS_MERGEDUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta) , path(input)
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
    def args      = task.ext.args  ?: ''
    def args2     = task.ext.args2 ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def extension = args2.contains("--output-fmt sam") ? "sam" :
                    args2.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    samtools merge \\
        ${args} \\
        - \\
        ${input} |\\
    samtools markdup \\
        -T ${prefix} \\
        -f ${prefix}.metrics \\
        --threads ${task.cpus} \\
        $reference \\
        $args2 \\
        - \\
        ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args  ?: ''
    def args2     = task.ext.args2 ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt sam") ? "sam" :
                    args2.contains("--output-fmt cram") ? "cram" :
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
