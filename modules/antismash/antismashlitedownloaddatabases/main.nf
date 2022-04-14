process ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::antismash-lite=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.0.1--pyhdfd78af_1' :
        'quay.io/biocontainers/antismash-lite:6.0.1--pyhdfd78af_1' }"
    containerOptions {
        workflow.containerEngine == 'singularity' ?
        "-B $database_css:/usr/local/lib/python3.8/site-packages/antismash/outputs/html/css,$database_detection:/usr/local/lib/python3.8/site-packages/antismash/detection,$database_modules:/usr/local/lib/python3.8/site-packages/antismash/modules" :
        workflow.containerEngine == 'docker' ?
        "-v $database_css:/usr/local/lib/python3.8/site-packages/antismash/outputs/html/css -v $database_detection:/usr/local/lib/python3.8/site-packages/antismash/detection -v $database_modules:/usr/local/lib/python3.8/site-packages/antismash/modules" :
        ''
        }
    input:
    path database_css
    path database_detection
    path database_modules

    output:
    path("antismash_db") , emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download-antismash-databases \\
        --database-dir antismash_db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
