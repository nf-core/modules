process AGAT_SPADDINTRONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/033434db0bd6ba28660401e1059286f36641fd8ce55faa11973fe5eaf312adcd/data' :
        'community.wave.seqera.io/library/agat:1.5.1--ae3cd948ce5e9795' }"

    input:
    tuple val(meta), path(gff)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/^v//'"), emit: versions_agat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ""
    output = "${prefix}.intron.gff"
    """
    agat_sp_add_introns.pl \\
        --gff ${gff} \\
        ${config_param} \\
        --out ${output} \\
        ${args}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.intron.gff"
    """
    touch ${output}

    """
}
