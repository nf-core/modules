process VIENNARNA_RNAFOLD {
    tag '$rna_fastq'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/viennarna:2.6.4--py310pl5321h6cc9453_1':
        'biocontainers/viennarna:2.6.4--py310pl5321h6cc9453_1' }"

    input:
    path fasta

    output:
    path "*.fold"         , emit: rnafold_txt
    path "*.ps"           , emit: rnafold_ps
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNAfold \\
        ${args} \\
        --jobs=${task.cpus} \\
        --infile=$fasta \\
        --outfile=${fasta}.fold

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNAfold: \$( RNAfold --version |& sed '1!d ; s/RNAfold //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${fasta}.fold
    touch ${fasta}.ps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNAfold: \$( RNAfold --version |& sed '1!d ; s/RNAfold //')
    END_VERSIONS
    """
}
