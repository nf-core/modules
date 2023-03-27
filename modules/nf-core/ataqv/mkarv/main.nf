process ATAQV_MKARV {
    label 'process_medium'

    conda "bioconda::ataqv=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ataqv:1.3.0--py39hccc85d7_2':
        'quay.io/biocontainers/ataqv:1.3.0--py39hccc85d7_2' }"

    input:
    path "jsons/*"

    output:
    path "html"        , emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkarv \\
        $args \\
        --concurrency $task.cpus \\
        --force \\
        ./html/ \\
        jsons/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        # mkarv: \$( mkarv --version ) # Use this when version string has been fixed
        ataqv: \$( ataqv --version )
    END_VERSIONS
    """
}
