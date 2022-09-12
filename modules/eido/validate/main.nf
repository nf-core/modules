process EIDO_VALIDATE {
    tag '$samplesheet'
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::eido=0.1.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/eido/0.1.9_cv1/eido_0.1.9_cv1.sif' :
        'biocontainers/eido:0.1.9_cv1' }"

    input:
    path samplesheet
    path schema

    output:
    path "versions.yml"           , emit: versions
    path "conversion.out"         , emit: conversion

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ $samplesheet = "https://" || $samplesheet = "http://" ]]; then
        SAMPLESHEET_PATH=$samplesheet
    else
        SAMPLESHEET_PATH=\$(readlink -f $samplesheet)
    fi

    eido validate $args --st-index sample \$SAMPLESHEET_PATH -s $schema -e > conversion.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //' ))
    END_VERSIONS
    """
}
