process TD2_PREDICT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ecc5fb9460ed4255919340ea296994a3491c91515e9583c842ed49599e22a9ce/data':
        "community.wave.seqera.io/library/td2:1.1.0--76046a413a4219c1"}"

    input:
    tuple val(meta), path(fasta), path(orfs_dir, stageAs: 'orfs')

    output:
    tuple val(meta), path("*.TD2.bed") , emit: bed
    tuple val(meta), path("*.TD2.cds") , emit: cds
    tuple val(meta), path("*.TD2.pep") , emit: pep
    tuple val(meta), path("*.TD2.gff3"), emit: gff3
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('TD2.Predict'), val("1.1.0"), emit: versions_td2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    TD2.Predict \\
        -t ${fasta} \\
        -O ${orfs_dir} \\
        --all-good \\
        --complete-orfs-only \\
        --verbose \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.TD2.bed
    touch ${prefix}.TD2.cds
    touch ${prefix}.TD2.pep
    touch ${prefix}.TD2.gff3
    """
}
