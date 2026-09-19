process VCF2ZARR_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e0bf4a8a21faa7319626812bc557404bb37b440df1af2bbc89a80771aca1f94/data':
        'community.wave.seqera.io/library/bio2zarr:0.1.7--a742d2d9b8ee4347' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcz"), emit: vcz
    tuple val("${task.process}"), val('vcf2zarr'), eval('vcf2zarr --version |& sed -n "s/.* //p"'), emit: versions_vcf2zarr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/bio2zarr/vcf2zarrconvert

Reason:
The module was incorrectly named. The installed package is bio2zarr, not vcf2zarr.
Per nf-core naming conventions, the module has been moved to bio2zarr/vcf2zarrconvert.
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
