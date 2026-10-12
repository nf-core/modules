process VCF2ZARR_EXPLODE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/7713d869c8c8259c10701c95fc105bad8bcfbd6735de6941a47c9e6e26e9bb2f/data':
        'community.wave.seqera.io/library/bio2zarr:0.1.8--c2c92dd3f64fb0f9' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.icf"), emit: icf
    tuple val("${task.process}"), val('vcf2zarr'), eval('vcf2zarr --version |& sed -n "s/.* //p"'), emit: versions_vcf2zarr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/bio2zarr/vcf2zarrexplode

Reason:
The module was incorrectly named. The installed package is bio2zarr, not vcf2zarr.
Per nf-core naming conventions, the module has been moved to bio2zarr/vcf2zarrexplode.
"""
    assert false: deprecation_message

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/bio2zarr/vcf2zarrconvert

Reason:
The module was incorrectly named. The installed package is bio2zarr, not vcf2zarr.
Per nf-core naming conventions, the module has been moved to bio2zarr/vcf2zarrconvert.
"""
    assert false: deprecation_message
}
