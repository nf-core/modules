process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.7--pyh7cba7a3_0' :
        'biocontainers/bwameth:0.2.7--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("BwamethIndex"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir BwamethIndex
    mv $fasta BwamethIndex/

    bwameth.py index BwamethIndex/$fasta

    rm BwamethIndex/$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir BwamethIndex
    touch BwamethIndex/$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """
}
