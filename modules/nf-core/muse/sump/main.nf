process MUSE_SUMP {
    tag "$meta.id"
    label 'process_medium'

    // TODO Update when maintainer publishes conda package and container
    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'docker.io/library/muse:2.1' }"

    container "docker.io/famkebaeuerle/muse:2.1"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MUSE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(muse_call_txt)
    tuple val(meta2), path(ref_vcf), path(ref_vcf_tbi)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // hands -G for WGS data and -E for WES data
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    /MuSE/MuSE \\
        sump \\
        $args \\
        -I $muse_call_txt \\
        -O ${prefix}.vcf \\
        -n $task.cpus    \\
        -D $ref_vcf

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
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muse: ${VERSION}
    END_VERSIONS
    """
}
