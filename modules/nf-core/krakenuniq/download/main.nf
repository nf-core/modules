process KRAKENUNIQ_DOWNLOAD {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2':
        'biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2' }"

    input:
    val pattern

    output:
    path "${pattern}/"  , emit: output
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    krakenuniq-download \\
        --threads ${task.cpus} \\
        -o ${pattern}/ \\
        ${pattern} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """
}
