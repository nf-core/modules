
process METABULI_BUILD {
    tag 'build'
    label 'process_medium'

    conda "bioconda::metabuli=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.5--pl5321h6a68c12_1':
        'biocontainers/metabuli:1.0.5--pl5321h6a68c12_1' }"

    input:
    path(db)

    output:
    path "metabuli_db.tar.gz", emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """

    ls $db/library/* > lib.txt
    metabuli \\
        build \\
        --threads $task.cpus \\
        $db \\
        lib.txt \\
        $acc2taxid \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli | grep Version | sed 's/^metabuli Version: //';))
    END_VERSIONS
    """

    stub:
    """
    touch metabuli_db.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli | grep Version | sed 's/^metabuli Version: //';))
    END_VERSIONS
    """
}
