process AGAT_SPMERGEANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/033434db0bd6ba28660401e1059286f36641fd8ce55faa11973fe5eaf312adcd/data' :
        'community.wave.seqera.io/library/agat:1.5.1--ae3cd948ce5e9795' }"

    input:
    tuple val(meta), path(gffs)
    path config

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    def file_names   = "${gffs}".split(' ')
    def gff_param    = file_names.collect { gff_file -> "--gff ${gff_file}" }.join(' ')
    if ( file_names.contains ( "${prefix}.gff" ) ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_merge_annotations.pl \\
        ${gff_param} \\
        ${config_param} \\
        ${args} \\
        --output ${prefix}.gff
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_names = "${gffs}".split(' ')
    if ( file_names.contains ( "${prefix}.gff" ) ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.gff
    """
}
