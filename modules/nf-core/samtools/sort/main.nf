process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)
    val index_format

    output:
    tuple val(meta), path("${prefix}.bam"),                 emit: bam,  optional: true
    tuple val(meta), path("${prefix}.cram"),                emit: cram, optional: true
    tuple val(meta), path("${prefix}.sam"),                 emit: sam,  optional: true
    tuple val(meta), path("${prefix}.${extension}.crai"),   emit: crai, optional: true
    tuple val(meta), path("${prefix}.${extension}.csi"),    emit: csi,  optional: true
    tuple val(meta), path("${prefix}.${extension}.bai"),    emit: bai,  optional: true
    path  "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    reference = fasta ? "--reference ${fasta}" : ""
    output_file = index_format ? "${prefix}.${extension}##idx##${prefix}.${extension}.${index_format} --write-index" : "${prefix}.${extension}"
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (extension == "sam") {
            error "Indexing not compatible with SAM output"
        }
    }
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${output_file} \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt cram") ? "cram" :
                "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (extension == "sam") {
            error "Indexing not compatible with SAM output"
        }
    }
    index = index_format ? "touch ${prefix}.${extension}.${index_format}" : ""

    """
    touch ${prefix}.${extension}
    ${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
