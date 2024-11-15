process MUSE_SUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/muse:2.1.2--f6ec9e78771509ff':
        'community.wave.seqera.io/library/muse:2.1.2--e8279641c6ef8c63' }"

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
    def VERSION = '2.1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cp $ref_vcf_tbi /tmp/
    rm $ref_vcf_tbi
    cp /tmp/$ref_vcf_tbi .
    rm /tmp/$ref_vcf_tbi

    MuSE \\
        sump \\
        $args \\
        -I $muse_call_txt \\
        -O ${prefix}.vcf \\
        -n $task.cpus    \\
        -D $ref_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: ${VERSION}
    END_VERSIONS
    """
}
