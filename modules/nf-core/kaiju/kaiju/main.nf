// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KAIJU_KAIJU {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/40cd79a6bb3b5feabe877bef8498d63a3b91fed44b30556af4c7e6554c4310f2/data'
        : 'community.wave.seqera.io/library/kaiju:1.10.3--8b4a2069cc29ac7a'}"

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

    def db_list = db instanceof List ? db : db.listDirectory()

    if (!db_list.find { db_files -> db_files.name.endsWith('names.dmp') } || !db_list.find { db_files -> db_files.name.endsWith('nodes.dmp') }) {
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
