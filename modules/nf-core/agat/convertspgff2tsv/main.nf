process AGAT_CONVERTSPGFF2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdf86a4db56d388b2cad1599b640e1e33ab964dcdb74a8147981d48a429bc7d5/data'
:         'community.wave.seqera.io/library/agat:1.6.1--1f8d94d325805bb4' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_convert_sp_gff2tsv.pl \\
        --gff ${gff} \\
        --output ${prefix}.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
