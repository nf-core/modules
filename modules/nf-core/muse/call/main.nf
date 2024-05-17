process MUSE_CALL {
    tag "$meta.id"
    label 'process_high'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'docker.io/library/muse:2.1' }"

    container "muse:2.1"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MUSE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.MuSE.txt"), emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    MuSE/MuSE call \\
        $args \\
        -f $reference \\
        -O ${prefix}  \\
        -n $task.cpus \\
        $tumor_bam    \\
        $normal_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muse: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.MuSE.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muse: ${VERSION}
    END_VERSIONS
    """
}
