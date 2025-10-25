process AGAT_SPFILTERBYORFSIZE {
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
    tuple val(meta), path("*.passed.gff"), emit: passed_gff
    tuple val(meta), path("*.failed.gff"), emit: failed_gff
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def config_arg  = config ? "-c ${config}" : ''
    if( "${gxf}" in [ "${prefix}.passed.gff", "${prefix}.failed.gff" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_filter_by_ORF_size.pl \\
        -g ${gxf} \\
        ${args} \\
        ${config_arg} \\
        -o ${prefix}

    mv \\
        ${prefix}_NOT* \\
        "${prefix}.failed.gff"

    mv \\
        ${prefix}_* \\
        "${prefix}.passed.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_flag_short_introns.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "${gxf}" in [ "${prefix}.passed.gff", "${prefix}.failed.gff" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.passed.gff
    touch ${prefix}.failed.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_flag_short_introns.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
