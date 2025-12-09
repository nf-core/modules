process TD2_LONGORFS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4155bf3b720e1e32d0615a947696fb0287ee4e8cdbeb4ec784dd4da7bb5b2e86/data':
        "community.wave.seqera.io/library/td2:1.0.6--ea3e5ac09443b677"}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/longest_orfs.{cds,gff3,pep}"), emit: orfs
    tuple val("${task.process}"), val('TD2.LongOrfs'), eval("echo ${VERSION}"), emit: versions_td2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    VERSION = '1.0.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    TD2.LongOrfs \\
    -t ${fasta} \\
    -O ${prefix} \\
    --threads ${task.cpus} \\
    ${args}
    """

    stub:
    VERSION = 'v1.0.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch ${prefix}/longest_orfs.cds
    touch ${prefix}/longest_orfs.gff3
    touch ${prefix}/longest_orfs.pep
    """
}
