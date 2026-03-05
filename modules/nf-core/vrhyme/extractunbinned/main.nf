process VRHYME_EXTRACTUNBINNED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vrhyme:1.1.0--pyhdfd78af_1':
        'biocontainers/vrhyme:1.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(membership)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fasta") , emit: unbinned_sequences
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_unbinned_sequences"
    def fasta_input = fasta.toString().replaceAll(/\.gz$/, '')
    def gunzip      = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${fasta_input}" : ""
    def cleanup     = fasta.getExtension() == "gz" ? "rm ${fasta_input}" : ""
    if( "$fasta_input" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    ${gunzip}

    extract_unbinned_sequences.py \\
        -i $membership \\
        -f $fasta_input \\
        -o ${prefix}.fasta \\
        $args

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_unbinned_sequences"
    def fasta_input = fasta.toString() - ~/\.gz$/
    if( "$fasta_input" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vrhyme: \$(echo \$(vRhyme --version 2>&1) | sed 's/^.*vRhyme v//; s/Using.*\$//' )
    END_VERSIONS
    """
}
