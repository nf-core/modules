process SENTIEON_APPLYVARCAL {
    tag "${meta.id}"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3aeac96f5787ca7dc1b8ed6679094246a6089ff8faf8f2d90d96fdace514229c/data'
        : 'community.wave.seqera.io/library/sentieon_gzip:b8170ec7a010f305'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(recal), path(recal_index), path(tranches)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi"),    emit: tbi
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_applyvarcal"
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver \\
        -r ${fasta}  \\
        -t ${task.cpus} \\
        ${args} \\
        --algo ApplyVarCal \\
        -v ${vcf} \\
        --recal ${recal} \\
        --tranches_file ${tranches} \\
        ${args2} \\
        ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_applyvarcal"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
