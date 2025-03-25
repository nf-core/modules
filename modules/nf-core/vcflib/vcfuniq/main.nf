process VCFLIB_VCFUNIQ {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/61/61442a6401d3dad42fc9645a00a4575420d306d345a9e9d694d031cf1b3f383f/data':
        'community.wave.seqera.io/library/vcflib:1.0.12--2281750e7717b014' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vcfuniq \\
        $vcf \\
        | bgzip -c $args > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: $VERSION
    END_VERSIONS
    """
}
