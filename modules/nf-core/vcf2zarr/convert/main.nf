process VCF2ZARR_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebfe707031ebecc7c4f597bb0f310465a6493ba44a16af56ab3d3872ee7492d2/data':
        'community.wave.seqera.io/library/bio2zarr:0.1.6--8ee007f51aad5560' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcz"), emit: vcz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2zarr \\
        convert \\
        $args \\
        --worker-processes $task.cpus \\
        $vcf \\
        ${prefix}.vcz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2zarr: \$(vcf2zarr --version |& sed -n 's/.* //p')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir ${prefix}.vcz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2zarr: \$(vcf2zarr --version |& sed -n 's/.* //p')
    END_VERSIONS
    """
}
