process MOTUS_PREPLONG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0':
        'biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.gz"), emit: out
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def refdb = db ? "-db ${db}" : ""

    """
    motus \\
        prep_long \\
        $args \\
        -i $reads \\
        $refdb \\
        -t $task.cpus \\
        -o ${prefix}.gz \\
        2> >(tee ${prefix}.log >&2)


    ## mOTUs version number is not available from command line.
    ## mOTUs save the version number in index database folder.
    ## mOTUs will check the database version is same version as exec version.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$(grep motus db_mOTU/db_mOTU_versions | sed 's/motus\\t//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$(grep motus db_mOTU/db_mOTU_versions | sed 's/motus\\t//g')
    END_VERSIONS
    """
}
