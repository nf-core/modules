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
    tuple val(meta), path("${prefix}/*.TD2.{bed,cds,gff3,pep}"), emit: predictions
    path("versions.yml")                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/

    TD2.Predict \\
        -t ${fasta} \\
        -O ${orfs_dir} \\
        ${args}

    mv *.TD2.{bed,cds,gff3,pep} ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        td2: \$(td2 v1.0.6)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch ${prefix}/fakefile.TD2.bed
    touch ${prefix}/fakefile.TD2.cds
    touch ${prefix}/fakefile.TD2.gff3
    touch ${prefix}/fakefile.TD2.pep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        td2: \$(td2 v1.0.6)
    END_VERSIONS
    """
}
