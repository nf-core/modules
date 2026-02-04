process TD2_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4155bf3b720e1e32d0615a947696fb0287ee4e8cdbeb4ec784dd4da7bb5b2e86/data':
        "community.wave.seqera.io/library/td2:1.0.6--ea3e5ac09443b677"}"

    input:
    tuple val(meta), path(fasta), path(orfs_dir, stageAs: 'orfs')

    output:
    tuple val(meta), path("${prefix}/*.TD2.{bed,cds,gff3,pep}")        , emit: predictions
    tuple val("${task.process}"), val('TD2.Predict'), eval("echo ${VERSION}"), emit: versions_td2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    VERSION = '1.0.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/

    TD2.Predict \\
        -t ${fasta} \\
        -O ${orfs_dir} \\
        ${args}

    mv *.TD2.{bed,cds,gff3,pep} ${prefix}/
    """

    stub:
    VERSION = 'v1.0.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch ${prefix}/${fasta}.TD2.bed
    touch ${prefix}/${fasta}.TD2.cds
    touch ${prefix}/${fasta}.TD2.gff3
    touch ${prefix}/${fasta}.TD2.pep
    """
}
