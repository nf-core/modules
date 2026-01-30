process SENTIEON_READWRITER {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}"),                             emit: output
    tuple val(meta), path("${prefix}.${index}"),                    emit: index
    tuple val(meta), path("${prefix}"), path("${prefix}.${index}"), emit: output_index
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:


    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def input_str = input.sort {in -> in.getName() }.collect {in ->  "-i ${in}" }.join(' ')
    def reference = fasta ? "-r ${fasta}" : ''

    // bam -> bam: prefix = "<filename>.bam"
    // bam -> cram: prefix = "<filename>.cram"
    // cram -> cram: prefix = "<filename>.cram"
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    index = prefix.tokenize('.')[-1] == "bam" ? "bai" : "crai"

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon \\
        driver \\
        -t ${task.cpus} \\
        ${reference} \\
        ${args} \\
        ${input_str} \\
        --algo ReadWriter \\
        ${args2} \\
        ${prefix}

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.cram"
    index = prefix.tokenize('.')[-1] == "bam" ? "bai" : "crai"
    """

    touch ${prefix}
    touch ${prefix}.${index}

    """
}
