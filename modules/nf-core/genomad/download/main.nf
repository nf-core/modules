process GENOMAD_DOWNLOAD {
    label 'process_single'

    conda "bioconda::genomad=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.5.2--pyhdfd78af_0':
        'biocontainers/genomad:1.5.2--pyhdfd78af_0' }"

    output:
    path "genomad_db/"  , emit: genomad_db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    genomad \\
        download-database .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
