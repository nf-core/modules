process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b542bbe1f99afd494ef07423ea8b52f2b8a081b85f92db2726c283c78da3cf0/data':
        'community.wave.seqera.io/library/biscuit_samtools:84373c8a97fa63b8' }"

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
