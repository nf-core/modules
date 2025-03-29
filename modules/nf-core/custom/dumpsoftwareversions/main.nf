def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This module is no longer recommended for use, as it is replaced by the function softwareVersionsToYAML
in the utils_nfcore_pipeline subworkflow that is included in the nf-core template.

"""
process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.27--pyhdfd78af_0' :
        'biocontainers/multiqc:1.27--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert true: deprecation_message
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
