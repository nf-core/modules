
process METABULI_BUILD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.5--pl5321h6a68c12_1':
        'biocontainers/metabuli:1.0.5--pl5321h6a68c12_1' }"

    input:
    tuple val(meta), path(db)
    path accession2taxid, stageAs: 'taxonomy/acc2taxid'

    output:
    tuple val(meta), path("$prefix"), emit: db
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    acc2taxid = accession2taxid ? "${accession2taxid}" : "${db}/taxonomy/acc2taxid"
    """
    find ${db}/library -type f -name '*.fna' > library-files.txt

    metabuli \\
        build \\
        ${db} \\
        library-files.txt \\
        ${acc2taxid} \\
        --db-name ${prefix} \\
        --threads ${task.cpus} \\
        $args

    if [[ \$(basename ${db}) != "${prefix}" ]]; then
        mkdir -p ${prefix}
        mv ${db}/* ${prefix}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """
}
