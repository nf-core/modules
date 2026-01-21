process SENTIEON_TNFILTER {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1dfe59ef66d7326b43db9ab1f39ce6220b358a311078c949a208f9c9815d4e/data'
        : 'community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(stats), path(contamination), path(segments), path(orientation_priors)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // options for the driver
    def args2 = task.ext.args2 ?: '' // options for --algo TNfilter
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    def contamination_command = contamination ? " --contamination ${contamination} " : ''
    def segments_command = segments ? segments.collect { "--tumor_segments ${it}" }.join(' ') : ''
    def orientation_priors_command = orientation_priors ? orientation_priors.collect { "--orientation_priors ${it}" }.join(' ') : ''
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver -r ${fasta} \\
    ${args} \\
    --algo TNfilter \\
    ${args2} \\
    -v ${vcf} \\
    ${contamination_command} \\
    ${segments_command} \\
    ${orientation_priors_command} \\
    ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    """
}
