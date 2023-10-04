process BEDTOOLS_JACCARD {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(input_a), path(input_b)
    tuple val(meta2), path(genome_file)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = genome_file ? "-g ${genome_file}" : ""
    """
    bedtools jaccard \\
        -a ${input_a} \\
        -b ${input_b} \\
        ${genome} \\
        ${args} \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
