process KRAKENUNIQ_DOWNLOAD {
    label 'process_low'

    conda "bioconda::krakenuniq=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.2--pl5321h19e8d03_0':
        'biocontainers/krakenuniq:1.0.2--pl5321h19e8d03_0' }"

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
