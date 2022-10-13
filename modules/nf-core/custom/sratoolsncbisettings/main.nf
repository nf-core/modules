process CUSTOM_SRATOOLSNCBISETTINGS {
    tag 'ncbi-settings'
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    def container_image = "/sra-tools:2.11.0--pl5321ha49a11a_3"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    output:
    path('*.mkfg')     , emit: ncbi_settings
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    template 'detect_ncbi_settings.sh'
}
