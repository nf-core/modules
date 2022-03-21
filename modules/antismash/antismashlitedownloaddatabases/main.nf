process ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES {
    tag '$database_dir'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::antismash-lite=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash-lite:6.0.1--pyhdfd78af_0' }"

    input:
    val(database_dir)

    output:
    path(database_dir) , emit: database

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def databases = database_dir ? "--database-dir ${database_dir}" : ''
    """
    download-antismash-databases \\
        $databases \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
