process PARAGRAPH_VCF2PARAGRAPH {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paragraph:2.3--h21f15d8_1':
        'biocontainers/paragraph:2.3--h21f15d8_1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.json.gz")  , emit: graph
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vcf2paragraph.py \\
        ${args} \\
        --reference-sequence ${fasta} \\
        ${vcf} \\
        ${prefix}.json

    bgzip --threads ${task.cpus} ${args2} ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo | gzip > ${prefix}.json.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
