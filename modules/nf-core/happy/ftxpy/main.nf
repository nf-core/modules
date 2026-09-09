process HAPPY_FTXPY {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.15--py27hcb73b3d_0':
        'quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0' }"

    input:
    tuple val(meta), path(vcf), path(regions_bed), path(targets_bed), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.csv"), emit: features
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('happy'), val('0.3.15'), topic: versions, emit: versions_happy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-R ${regions_bed}" : ""
    def targets = targets_bed ? "-T ${targets_bed}" : ""
    def bams = bam ? "--bam ${bam}" : ""

    """
    ftx.py \\
        -o ${prefix} \\
        ${regions} \\
        ${targets} \\
        ${bams} \\
        --reference ${fasta} \\
        ${args} \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
