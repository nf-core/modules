process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda "bioconda::sra-tools=3.0.8 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:2f4a4c900edd6801ff0068c2b3048b4459d119eb-0' :
        'biocontainers/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:2f4a4c900edd6801ff0068c2b3048b4459d119eb-0' }"

    input:
    tuple val(meta), val(id)
    path ncbi_settings
    path certificate

    output:
    tuple val(meta), path("${id}{,.sralite}"), emit: sra
    path 'versions.yml'                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = task.ext.args ?: ''
    args += certificate ? " --perm ${certificate}" : ''
    args2 = task.ext.args2 ?: '5 1 100'  // <num retries> <base delay in seconds> <max delay in seconds>
    template 'retry_with_backoff.sh'
}
