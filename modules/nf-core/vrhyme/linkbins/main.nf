process VRHYME_LINKBINS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vrhyme=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1':
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*_linked_bins.fasta")        , emit: linked_bins
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    link_bin_sequences.py \\
        -i $bins \\
        -o vRhyme_linked_bins \\
        -e fasta \\
        $args

    cat vRhyme_linked_bins/*.fasta > ${prefix}_linked_bins.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_linked_bins.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """
}
