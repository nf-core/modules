process VCF2DB {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/3013992b36b50c203acfd01b000d37f3753aee640238f6dd39d5e47f58e54d98/data':
        'community.wave.seqera.io/library/python_python-snappy_snappy_cyvcf2_vcf2db:9c1d7f361187f21a' }"

    input:
    tuple val(meta), path(vcf), path(ped)

    output:
    tuple val(meta), path("*.db") , emit: db
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('vcf2db'), val("2020.02.24"), emit: versions_vcf2db, topic: versions
    tuple val("${task.process}"), val('python-snappy'), val("0.5.4"), emit: versions_python_snappy, topic: versions
    tuple val("${task.process}"), val('snappy'), val("1.1.8"), emit: versions_snappy, topic: versions
    tuple val("${task.process}"), val('cyvcf2'), eval("python -c 'import cyvcf2; print(cyvcf2.__version__)'"), emit: versions_cyvcf2, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | awk '{print \$2}'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2db.py \\
        $vcf \\
        $ped \\
        ${prefix}.db \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.db
    """
}
