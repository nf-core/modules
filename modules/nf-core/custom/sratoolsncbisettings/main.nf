process CUSTOM_SRATOOLSNCBISETTINGS {
    tag 'ncbi-settings'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.2.1--h4304569_1' :
        'biocontainers/sra-tools:3.2.1--h4304569_1' }"

    input:
    val ids

    output:
    path('*.mkfg'), emit: ncbi_settings
    tuple val("${task.process}"), val('sratools'), eval("prefetch --version 2>&1 | grep -Eo '[0-9.]+'"), topic: versions, emit: versions_sratools

    when:
    task.ext.when == null || task.ext.when

    script:
    config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"

    """
    echo ${config}
    """

    template "detect_ncbi_settings.sh"

    stub:
    """
    touch user-settings.mkfg
    """
}
