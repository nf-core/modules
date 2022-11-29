process CHECKV_DOWNLOADDATABASE {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::checkv=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
    path fasta
    path db

    output:
    path "checkv-db-*"           , emit: checkv_db
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def checkv_db = db ?: ''
    def update_sequence = fasta ?: ''

    // determine which command needs to be executed
    if (checkv_db != '' & update_sequence != '') {
        method = "checkv update_database --threads $task.cpus"
    }else if (checkv_db != '' & update_sequence == '') {
        method = 'ln -s'
    }else{
        method = 'checkv download_database'
    }

    """
    $method \\
        $args \\
        $db \\
        ./  \\
        $fasta \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

}
