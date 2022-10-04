process CUSTOM_SRATOOLSNCBISETTINGS {
    tag 'ncbi-settings'
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5321ha49a11a_3' :
        'quay.io/biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3' }"

    output:
    path('*.mkfg')     , emit: ncbi_settings
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    template 'detect_ncbi_settings.sh'
}
