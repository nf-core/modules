// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MIRTRACE_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirtrace:1.0.1--0':
        'quay.io/biocontainers/mirtrace:1.0.1--0' }"

    input:
    tuple val(meta), path(reads), path(mirtrace_config)
    val(mirtrace_species)

    output:
    tuple val(meta), path ("*.html")                                                , emit: html
    tuple val(meta), path ("*.json")                                                , emit: json
    tuple val(meta), path ("*.tsv")                                                 , emit: tsv
    tuple val(meta), path ("qc_passed_reads.all.collapsed/*.{fa,fasta}")            , emit: all_fa
    tuple val(meta), path ("qc_passed_reads.rnatype_unknown.collapsed/*.{fa,fasta}"), emit: rnatype_unknown_fa
    tuple val("${task.process}"), val('mirtrace'), eval("mirtrace -v"), emit: versions_mirtrace, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mirtrace_mode = mirtrace_config ? "--config ${mirtrace_config}": "${reads}"

    """
    mirtrace qc  \\
        --species ${mirtrace_species} \\
        --write-fasta \\
        --output-dir . \\
        --force \\
        ${mirtrace_mode}
    """

    stub:
    // Mirror the real output names: miRTrace names the collapsed FASTA files after the input
    // read files (without the FASTQ extension) and always writes 7 statistics TSV files
    def fasta_names = [reads].flatten().collect { read -> read.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '') + '.fasta' }
    """
    touch mirtrace-report.html
    touch mirtrace-results.json
    for stats in contamination_basic contamination_detailed length mirna-complexity phred qcstatus rnatype; do
        touch mirtrace-stats-\${stats}.tsv
    done

    mkdir -p qc_passed_reads.all.collapsed
    mkdir -p qc_passed_reads.rnatype_unknown.collapsed

    for fasta in ${fasta_names.join(' ')}; do
        touch qc_passed_reads.all.collapsed/\${fasta}
        touch qc_passed_reads.rnatype_unknown.collapsed/\${fasta}
    done
    """
}
