process AGAT_SPKEEPLONGESTISOFORM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/033434db0bd6ba28660401e1059286f36641fd8ce55faa11973fe5eaf312adcd/data' :
        'community.wave.seqera.io/library/agat:1.5.1--ae3cd948ce5e9795' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ""
    output           = "${prefix}.longest.gff"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff ${gxf} \\
        ${config_param} \\
        --out ${output} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output     = "${prefix}.longest.gff"
    """
    touch ${output}
    touch ${gxf}.agat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat --version)
    END_VERSIONS
    """
}
