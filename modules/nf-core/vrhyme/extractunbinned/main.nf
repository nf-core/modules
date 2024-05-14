process VRHYME_EXTRACTUNBINNED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vrhyme=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1':
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(membership)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_unbinned_sequences.fasta") , emit: unbinned_sequences
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_unbinned_sequences.py \\
        -i $membership \\
        -f $fasta \\
        -o ${prefix}_unbinned_sequences.fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unbinned_sequences.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """
}
