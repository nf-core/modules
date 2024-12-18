process CUSTOM_SRATOOLSNCBISETTINGS {
    tag 'ncbi-settings'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.0.8--h9f5acd7_0' :
        'biocontainers/sra-tools:3.0.8--h9f5acd7_0' }"

    input:
    val ids

    output:
    path("*.mkfg")     , emit: ncbi_settings
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    template 'detect_ncbi_settings.sh'
}
