process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda (params.enable_conda ? 'bioconda::multiqc=1.13' : null)
        def container_image = "/multiqc:1.13--pyhdfd78af_0"
                                                              container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
