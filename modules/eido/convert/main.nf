process EIDO_CONVERT {
    tag '$samplesheet'
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::eido=0.1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/eido/0.1.9_cv1/eido_0.1.9_cv1.sif' :
        'biocontainers/eido:0.1.9_cv1' }"

    input:
    path samplesheet

    output:
    path "versions.yml"           , emit: versions
    path "*.csv"                  , emit: samplesheet_converted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    eido \\
        convert \\
        $samplesheet \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' ))
    END_VERSIONS
    """
}
