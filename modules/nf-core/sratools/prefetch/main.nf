process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda "bioconda::sra-tools=2.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5321ha49a11a_3' :
        'biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3' }"

    input:
    tuple val(meta), val(id)
    path ncbi_settings
    path certificate

    output:
    tuple val(meta), path(id), emit: sra
    path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = task.ext.args ?: ''
    args += certificate ? " --perm ${certificate}" : ''
    args2 = task.ext.args2 ?: '5 1 100'  // <num retries> <base delay in seconds> <max delay in seconds>
    template 'retry_with_backoff.sh'
}
