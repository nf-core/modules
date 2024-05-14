process PURGEDUPS_GETSEQS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::purge_dups=1.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0':
        'biocontainers/purge_dups:1.2.6--h7132678_0' }"

    input:
    tuple val(meta), path(assembly), path(bed)

    output:
    tuple val(meta), path("*.hap.fa")   , emit: haplotigs
    tuple val(meta), path("*.purged.fa"), emit: purged
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_seqs \\
        $args \\
        -e $bed \\
        -p $prefix \\
        $assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """
}
