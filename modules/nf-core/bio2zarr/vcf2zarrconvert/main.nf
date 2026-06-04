process BIO2ZARR_VCF2ZARRCONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/7713d869c8c8259c10701c95fc105bad8bcfbd6735de6941a47c9e6e26e9bb2f/data':
        'community.wave.seqera.io/library/bio2zarr:0.1.8--c2c92dd3f64fb0f9' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcz"), emit: vcz
    tuple val("${task.process}"), val('vcf2zarr'), eval('vcf2zarr --version |& sed -n "s/.* //p"'), topic: versions  , emit: versions_vcf2zarr

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
