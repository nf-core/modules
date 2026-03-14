process VCF2ZARR_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2zarr \\
        convert \\
        ${args} \\
        --worker-processes ${task.cpus} \\
        ${vcf} \\
        ${prefix}.vcz
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir ${prefix}.vcz
    """
}
