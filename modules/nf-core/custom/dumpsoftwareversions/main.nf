// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/043d2f0d10cc53f742deb85a68f02fd5edb2c450a40eb6bf8ecfd40ee95b7a02/data'
        : 'community.wave.seqera.io/library/multiqc:1.27--5aad79bfb397e7db'}"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This module is no longer recommended for use, as it is replaced by the function softwareVersionsToYAML
in the utils_nfcore_pipeline subworkflow that is included in the nf-core template.

"""
    assert false: deprecation_message
    template 'dumpsoftwareversions.py'
}
