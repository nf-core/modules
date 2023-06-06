process PARAGRAPH_IDXDEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::paragraph=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paragraph:2.3--h8908b6f_0':
        'biocontainers/paragraph:2.3--h8908b6f_0' }"

    input:
    tuple val(meta) , path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.json") , emit: depth
    tuple val(meta), path("*.tsv")  , emit: binned_depth, optional:true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def type = input.extension
    def output_bins = type == "cram" ? "--output-bins ${prefix}.tsv" : ""
    """
    idxdepth \\
        --bam ${input} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --output ${prefix}.json \\
        ${output_bins} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def type = input.extension
    def output_bins = type == "cram" ? "touch ${prefix}.tsv" : ""
    """
    touch ${prefix}.json
    ${output_bins}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
    END_VERSIONS
    """
}
