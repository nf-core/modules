process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.1.0.20220707--he272189_1':
        'biocontainers/biscuit:1.1.0.20220707--he272189_1' }"

    input:
    path fasta, stageAs: "BiscuitIndex/*"

    output:
    path "BiscuitIndex/*.fa*", emit: index, includeInputs: true
    path "versions.yml"      , emit: versions

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
