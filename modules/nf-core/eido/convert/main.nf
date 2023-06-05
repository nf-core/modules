process EIDO_CONVERT {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::eido=0.1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/eido/0.1.9_cv1/eido_0.1.9_cv1.sif' :
        'docker.io/biocontainers/eido:0.1.9_cv1' }"

    input:
    path samplesheet
    val format
    path pep_input_base_dir

    output:
    path "versions.yml"           , emit: versions
    path "${prefix}.${format}"    , emit: samplesheet_converted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "samplesheet_converted"
    """
    eido \\
        convert \\
        -f $format \\
        $samplesheet \\
        $args \\
        -p samples=${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' )
    END_VERSIONS
    """
}
