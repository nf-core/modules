process SAMTOOLS_MARKDUP {
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
    tuple val(meta), path("*.bam"),  emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.sam"),  emit: sam,  optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools \\
        markdup \\
        $args \\
        ${reference} \\
        -@ $task.cpus \\
        -T $prefix \\
        $input \\
        ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}
    """
}
