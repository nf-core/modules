process CENTRIFUGER_QUANTIFICATION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuger:1.1.0--hf426362_0':
        'quay.io/biocontainers/centrifuger:1.1.0--hf426362_0' }"

    input:
    tuple val(meta), path(classification_file)
    tuple val(meta2), path(db)
    path taxonomy_nodes
    path taxonomy_names
    path size_table

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: report_file
    tuple val("${task.process}"), val('centrifuger'), eval("centrifuger -v 2>&1 | sed 's/Centrifuger v//'"), emit: versions_centrifuger,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // include -x option with index prrefix or use specified files
    def db_arg = ""
    if (db) {
        db_arg= " -x `find -L ${db} -name '*.1.cfr' -not -name '._*'  | sed 's/\\.1.cfr\$//'`"
        }
    else {
        def tax_arg = taxonomy_nodes ? "--taxonomy-tree ${taxonomy_nodes}" : ""
        def name_arg = taxonomy_names ? "--name-table ${taxonomy_names}" : ""
        def size_arg = size_table ? "--size-table ${size_table}" : ""
        db_arg = "${tax_arg} ${name_arg} ${size_arg}"
        }

    """
    centrifuger-quant \\
        ${db_arg} \\
        -c ${classification_file} \\
        ${args} > ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    #output
    echo "" > ${prefix}.tsv
    """
}
