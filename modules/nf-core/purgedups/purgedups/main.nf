process PURGEDUPS_PURGEDUPS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::purge_dups=1.2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0':
        'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0' }"

    input:
    tuple val(meta), path(basecov), path(cutoff), path(paf)

    output:
    tuple val(meta), path("*.dups.bed")      , emit: bed
    tuple val(meta), path("*.purge_dups.log"), emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    purge_dups \\
        $args \\
        -T $cutoff \\
        -c $basecov \\
        $paf > ${prefix}.dups.bed 2> ${prefix}.purge_dups.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """
}
