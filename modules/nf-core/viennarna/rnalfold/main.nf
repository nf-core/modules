process VIENNARNA_RNALFOLD {
    tag '$rna_fastq'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/viennarna:2.6.4--py310pl5321h6cc9453_1':
        'quay.io/biocontainers/viennarna:2.6.4--py310pl5321h6cc9453_1' }"

    input:
    path fasta

    output:
    path "*.lfold"                                                                       , emit: rnalfold_txt
    tuple val("${task.process}"), val('RNALfold'), eval("RNALfold --version 2>&1 | sed -n '1s/RNALfold //p'"), emit: versions_rnalfold, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RNALfold \\
        ${args} \\
        --infile=$fasta \\
        --outfile=${fasta.baseName}.lfold
    """

    stub:
    """
    touch ${fasta.baseName}.lfold
    touch ${fasta}.ps
    """
}
