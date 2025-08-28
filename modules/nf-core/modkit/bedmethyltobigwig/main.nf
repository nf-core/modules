process MODKIT_BEDMETHYLTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30413db983a1908edc096dcd8a31e7c26363b2010c5293e6bf59c0cef5d4a88a/data':
        'biocontainers/ont-modkit:0.5.0--hcdda2d0_1' }"

    input:
    tuple val(meta),  path(bedmethyl)
    tuple val(meta2), path(chromsizes)
    val modcodes

    output:
    tuple val(meta), path("*.bw"), emit: bw
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mods = modcodes instanceof List ? modcodes.join(',') : modcodes
    """
    modkit bedmethyl tobigwig \\
        $args \\
        --nthreads $task.cpus \\
        --sizes $chromsizes \\
        --mod-codes $mods \\
        $bedmethyl \\
        ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """
}
