process BLOBTK_CREATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/08833d1b91f41024e06e2cb5a982598531199c04e6544885d42ef2cb0480de18/data' :
        'community.wave.seqera.io/library/blobtk:0.8.0--2fe0d833a26e0cd9' }"

    input:
    tuple val(meta), path(fasta), path(full_table)  // [meta, fasta, full_table]

    output:
    tuple val(meta), path("${prefix}/"), emit: blobdir
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix              = task.ext.prefix ?: "${meta.id}"

    def full_table_args = full_table ? "--busco ${full_table}"  : ""

    """
    blobtk create \\
        --fasta ${fasta} \\
        ${full_table_args} \\
        --out ${prefix}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/meta.json
    touch ${prefix}/gc.json
    touch ${prefix}/identifiers.json
    touch ${prefix}/n.json
    touch ${prefix}/ncount.json
    touch ${prefix}/length.json
    """
}
