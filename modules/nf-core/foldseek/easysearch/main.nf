process FOLDSEEK_EASYSEARCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek:10.941cd33--h5021889_1':
        'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1' }"

    input:
    tuple val(meta) , path(pdb)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("${prefix}.m8"), emit: aln
    tuple val("${task.process}"), val('foldseek'), eval("foldseek --help |& sed -n 's/.*Version: //p'"), emit: versions_foldseek, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB_NAME="${db}/${meta2.id}"
    if [ ! -f "\${DB_NAME}.dbtype" ]; then
        DB_NAME=\$(find -L "${db}/" -maxdepth 1 -name '*.lookup' | sed 's/\\.lookup\$//')
        if [ "\$(printf '%s\\n' "\$DB_NAME" | grep -c .)" -ne 1 ]; then
            echo "ERROR: ${db}/${meta2.id} not found and expected exactly one .lookup file in ${db}" >&2
            exit 1
        fi
    fi

    foldseek \\
        easy-search \\
        ${pdb} \\
        "\$DB_NAME" \\
        ${prefix}.m8 \\
        tmp \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.m8
    """
}
