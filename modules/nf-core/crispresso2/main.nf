process CRISPRESSO2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crispresso2:2.3.4--py312hfcd9dac_0' :
        'quay.io/biocontainers/crispresso2:2.3.4--py312hfcd9dac_0' }"

    input:
    tuple val(meta), path(reads)
    val amplicon_sequences
    path amplicon_file

    output:
    tuple val(meta), path("CRISPResso_on_*")        , emit: results
    tuple val(meta), path("*.html")                 , emit: html
    tuple val(meta), path("CRISPResso_on_*/*.txt")  , emit: txt
    tuple val("${task.process}"), val('crispresso2'), eval("CRISPResso --version 2>&1 | sed 's/CRISPResso //'"), emit: versions_crispresso2, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle single-end vs paired-end reads
    def read_inputs = ""
    if (reads instanceof Path || reads.size() == 1) {
        read_inputs = "-r1 ${reads}"
    } else {
        read_inputs = "-r1 ${reads[0]} -r2 ${reads[1]}"
    }

    // Handle amplicon sequence vs amplicon file
    def amplicon_input = ""
    if (amplicon_file && amplicon_sequences) {
        error "Both amplicon_file and amplicon_sequences are provided. Please provide only one."
    } else if (amplicon_file) {
        amplicon_input = "-a \"\$(paste -sd, ${amplicon_file})\""
    } else if (amplicon_sequences) {
        amplicon_input = "-a ${amplicon_sequences}"
    } else if (task.ext.args && !task.ext.args.contains("--auto")) {
        error "Neither amplicon_file nor amplicon_sequences is provided and automatic amplicon detection is not enabled. Please provide one."
    }

    """
    export MPLCONFIGDIR=.matplotlib
    CRISPResso \\
        ${read_inputs} \\
        ${amplicon_input} \\
        --name ${prefix} \\
        --output_folder . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env bash
    set -e
    export MPLCONFIGDIR=.matplotlib

    mkdir -p CRISPResso_on_${prefix}
    touch CRISPResso_on_${prefix}/${prefix}.html
    touch CRISPResso_on_${prefix}/CRISPResso_report.html
    touch CRISPResso_on_${prefix}/CRISPResso_mapping_statistics.txt
    touch CRISPResso_on_${prefix}/quantification_of_editing_frequency.txt
    touch ${prefix}.html
    """
}
