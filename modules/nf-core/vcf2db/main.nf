process VCF2DB {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14d1257c98f789e23a888e2961673b5b98d89e4d03e6a3efba2b1134ed439f61/data':
        'community.wave.seqera.io/library/python_vcf2db:91f604106ada5cf2' }"

    input:
    tuple val(meta), path(vcf), path(ped)

    output:
    tuple val(meta), path("*.db") , emit: db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2020.02.24" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vcf2db.py \\
        $vcf \\
        $ped \\
        ${prefix}.db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2db: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2020.02.24" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2db: $VERSION
    END_VERSIONS
    """
}
