process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"),     emit: bam,  optional: true
    tuple val(meta), path("*.cram"),    emit: cram, optional: true
    tuple val(meta), path("*.crai"),    emit: crai, optional: true
    tuple val(meta), path("*.csi"),     emit: csi,  optional: true
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        --threads $task.cpus \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
