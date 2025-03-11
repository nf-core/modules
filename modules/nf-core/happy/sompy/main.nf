process HAPPY_SOMPY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.14--py27h5c5a3ab_0':
        'biocontainers/hap.py:0.3.14--py27h5c5a3ab_0' }"

    input:
    tuple val(meta), path(query_vcf), path(truth_vcf), path(regions_bed), path(targets_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(false_positives_bed)
    tuple val(meta5), path(ambiguous_beds)
    tuple val(meta6), path(bams)

    output:
    tuple val(meta), path('*.features.csv')           , emit: features, optional: true
    tuple val(meta), path('*.metrics.json')           , emit: metrics
    tuple val(meta), path('*.stats.csv')              , emit: stats
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-R ${regions_bed}" : ""
    def targets = targets_bed ? "-T ${targets_bed}" : ""
    def false_positives = false_positives_bed ? "--false-positives ${false_positives_bed}" : ""
    def ambiguous = ambiguous_beds ? "--ambiguous ${ambiguous_beds}" : ""
    def bams = bams ? "--bam ${bams}" : ""
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    som.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        ${args} \\
        --reference ${fasta} \\
        ${regions} \\
        ${targets} \\
        ${false_positives} \\
        ${ambiguous} \\
        ${bams} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.features.csv
    touch ${prefix}.metrics.json
    touch ${prefix}.stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """
}
