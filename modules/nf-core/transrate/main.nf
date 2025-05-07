process TRANSRATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/55/552b9db550aad1547ac928314812e788cc516f3ef6f6782ddc7dffa03c4d2145/data':
        'community.wave.seqera.io/library/snap-aligner_transrate:501869af4f81472b' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}"), emit: assembly
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    transrate \\
        $args \\
        ${fastq} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: \$(transrate --version 2>/dev/null | tail -n1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    touch ${prefix}/test.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: \$(transrate --version 2>/dev/null | tail -n1)
    END_VERSIONS
    """
}
