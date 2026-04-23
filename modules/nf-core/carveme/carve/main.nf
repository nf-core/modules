process CARVEME_CARVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/carveme:1.6.6--pyhdfd78af_1'
        : 'quay.io/biocontainers/carveme:1.6.6--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(fasta), path(universe_file), path(mediadb), path(soft), path(hard), path(reference)

    output:
    tuple val(meta), path("${prefix}.xml"), emit: model
    tuple val("${task.process}"), val('carveme'), eval("pip show carveme | sed -n 's/^Version: //p'"), topic: versions, emit: versions_carveme

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def universe_file_arg = universe_file ? "--universe-file ${universe_file}" : ''
    def mediadb_arg       = mediadb       ? "--mediadb ${mediadb}"             : ''
    def soft_arg          = soft          ? "--soft ${soft}"                   : ''
    def hard_arg          = hard          ? "--hard ${hard}"                   : ''
    def reference_arg     = reference     ? "--reference ${reference}"         : ''
    """
    carve \\
        ${fasta} \\
        --output ${prefix}.xml \\
        ${universe_file_arg} \\
        ${mediadb_arg} \\
        ${soft_arg} \\
        ${hard_arg} \\
        ${reference_arg} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.xml
    """
}
