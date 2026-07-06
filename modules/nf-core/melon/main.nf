process MELON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/melon:0.2.5--pyhdfd78af_0':
        'quay.io/biocontainers/melon:0.2.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path(database)
    path(k2_db)

    output:
    tuple val(meta), path("${prefix}/*.tsv") , emit: tsv_output
    tuple val(meta), path("${prefix}/*.json"), emit: json_output
    tuple val(meta), path("${prefix}.log")   , emit: log
    tuple val("${task.process}"), val("melon"), eval("melon -v"), emit: versions_melon, topic: versions
    tuple val("${task.process}"), val("diamond"), eval("diamond version | sed 's/.* version //'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def k2_db_arg = k2_db ? "--db-kraken ${k2_db}" : ''
    """
    melon \\
        ${reads} \\
        --db ${database} \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        ${k2_db_arg} \\
        ${args} \\
        2>| >(tee ${prefix}.log >&2)
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.json
    touch ${prefix}.log
    """
}
