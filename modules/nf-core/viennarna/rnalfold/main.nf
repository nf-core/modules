process VIENNARNA_RNALFOLD {
    tag '$rna_fastq'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/viennarna:2.6.4--py310pl5321h6cc9453_1':
        'biocontainers/viennarna:2.6.4--py310pl5321h6cc9453_1' }"

    input:
    path fasta

    output:
    path "*.lfold"        , emit: rnalfold_txt
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNALfold \\
        ${args} \\
        --infile=$fasta \\
        --outfile=${fasta.baseName}.lfold

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNALfold: \$( RNAfold --version |& sed '1!d ; s/RNALfold //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${fasta.baseName}.lfold
    touch ${fasta}.ps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNALfold: \$( RNAfold --version |& sed '1!d ; s/RNALfold //')
    END_VERSIONS
    """
}
