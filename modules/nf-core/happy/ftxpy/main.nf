process HAPPY_FTXPY {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.14--py27h5c5a3ab_0':
        'biocontainers/hap.py:0.3.14--py27h5c5a3ab_0' }"

    input:
    tuple val(meta), path(vcf), path(regions_bed), path(targets_bed), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.csv"), emit: features
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-R ${regions_bed}" : ""
    def targets = targets_bed ? "-T ${targets_bed}" : ""
    def bams = bam ? "--bam ${bam}" : ""
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ftx.py \\
        -o ${prefix} \\
        $regions \\
        $targets \\
        $bams \\
        --reference ${fasta} \\
        $args \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-R ${regions_bed}" : ""
    def targets = targets_bed ? "-T ${targets_bed}" : ""
    def bams = bam ? "--bam ${bam}" : ""
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """
}
