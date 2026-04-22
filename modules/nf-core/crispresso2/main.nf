process CRISPRESSO2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crispresso2:2.3.3--py39hff726c5_0' :
        'quay.io/biocontainers/crispresso2:2.3.3--py39hff726c5_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("CRISPResso_on_*")        , emit: results
    tuple val(meta), path("*.html")                 , emit: html
    tuple val(meta), path("CRISPResso_on_*/*.txt")  , emit: txt
    tuple val("${task.process}"), val('crispresso2'), eval("CRISPResso --version 2>&1 | sed 's/CRISPResso //'"), emit: versions_crispresso2, topic: versions

    when:
    task.ext.when == null || task.ext.when

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

    """
    export MPLCONFIGDIR=.matplotlib
    CRISPResso \\
        ${read_inputs} \\
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
