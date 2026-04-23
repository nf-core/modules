process KAIJU_KAIJU {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0'
        : 'quay.io/biocontainers/kaiju:1.10.0--h43eeafb_0'}"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path('*.tsv'), emit: results
    tuple val("${task.process}"), val('kaiju'), eval("kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'"), emit: versions_kaiju, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-i ${reads}" : "-i ${reads[0]} -j ${reads[1]}"

    if (!db.find { db_files -> db_files.name.endsWith('names.dmp') } || !db.find { db_files -> db_files.name.endsWith('nodes.dmp') }) {
        error('[KAIJU_KAIJU] Module error: Missing one of `nodes.dmp`, `names.dmp`. Check input.')
    }
    """
    dbnodes=`find -L ${db} -name "*nodes.dmp"`
    dbname=`find -L ${db} -name "*.fmi" -not -name "._*"`
    kaiju \\
        ${args} \\
        -z ${task.cpus} \\
        -t \$dbnodes \\
        -f \$dbname \\
        -o ${prefix}.tsv \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
