process MUSE_SUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0ebb574ef5eed2a6e034f1b2feea6c252d1ab0c8bc5135a669059aa1f4d2ca/data':
        'community.wave.seqera.io/library/muse:6637291dcbb0bdb8' }"

    input:
    tuple val(meta), path(muse_call_txt)
    tuple val(meta2), path(ref_vcf), path(ref_vcf_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '' // hands -G for WGS data and -E for WES data
    def args2  = task.ext.args2  ?: '' // args for gzip
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MuSE \\
        sump \\
        $args \\
        -I $muse_call_txt \\
        -n $task.cpus    \\
        -D $ref_vcf \\
        | gzip $args2 --stdout > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
    END_VERSIONS
    """
}
