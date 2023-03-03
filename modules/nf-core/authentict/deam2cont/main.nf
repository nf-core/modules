process AUTHENTICT_DEAM2CONT {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::authentict=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/authentict:1.0.0--py311h9f5acd7_1':
        'quay.io/biocontainers/authentict:1.0.0--py311h9f5acd7_1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(config)
    tuple val(meta3), path(positions)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def config_file = config ? "-c ${config}" : ""
    def positions_file = positions ? "-p ${positions}" : ""

    """
    samtools view $args $bam | AuthentiCT \\
        deam2cont \\
        $args2 \\
        $config_file \\
        $positions_file \\
        - \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        authentict: $VERSION
    END_VERSIONS
    """
}
