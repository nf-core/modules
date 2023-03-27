process VIBRANT_DOWNLOADDB {
    label 'process_medium'

    conda "bioconda::vibrant=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vibrant:1.2.1--hdfd78af_3':
        'quay.io/biocontainers/vibrant:1.2.1--hdfd78af_3' }"

    output:

    path "vibrant_db"           , emit: db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    download-db.sh \\
        vibrant_db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VIBRANT: \$(echo \$(VIBRANT_run.py --version 2>&1) | sed 's/^.*VIBRANT v//;')
    END_VERSIONS
    """
}
