process PARAGRAPH_GRMPY {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::paragraph=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paragraph:2.3--h8908b6f_0':
        'quay.io/biocontainers/paragraph:2.3--h8908b6f_0' }"

    input:
    tuple val(meta), path(variants), path(variants_index), path(cram), path(crai), path(manifest)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    multigrmpy.py \\
        --input ${variants} \\
        --manifest ${manifest} \\
        --output ${prefix} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
    END_VERSIONS
    """
}
