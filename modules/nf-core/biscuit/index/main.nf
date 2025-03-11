process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33a9ca30b4154f11253c8d91a75382065dcb8282ba99b74dbee59ed8faceabd7/data':
        'community.wave.seqera.io/library/biscuit:1.5.0.20240506--ca92d9d0a37b5fa8' }"

    input:
    tuple val(meta), path(fasta, name:"BiscuitIndex/")

    output:
    tuple val(meta), path("BiscuitIndex"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    biscuit \\
        index \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${fasta}.bis.amb
    touch ${fasta}.bis.ann
    touch ${fasta}.bis.pac
    touch ${fasta}.dau.bwt
    touch ${fasta}.dau.sa
    touch ${fasta}.par.bwt
    touch ${fasta}.par.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}
