process SENTIEON_BWAINDEX {
    tag "${fasta}"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "bwa/${task.ext.prefix}" : "bwa/${fasta.baseName}"
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}
    mkdir bwa

    sentieon \\
        bwa index \\
        ${args} \\
        -p ${prefix} \\
        ${fasta}
    """

    stub:
    """
    mkdir bwa

    touch bwa/genome.amb
    touch bwa/genome.ann
    touch bwa/genome.bwt
    touch bwa/genome.pac
    touch bwa/genome.sa
    """
}
