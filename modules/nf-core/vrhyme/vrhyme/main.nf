process VRHYME_VRHYME {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vrhyme=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1':
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("**/vRhyme_bin_*.fasta")                  , emit: bins
    tuple val(meta), path("**/vRhyme_best_bins.*.membership.tsv")   , emit: membership
    tuple val(meta), path("**/vRhyme_best_bins.*.summary.tsv")      , emit: summary
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vRhyme \\
        -i $fasta \\
        -b $reads \\
        -o vRhyme \\
        -t $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vrhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """
}
