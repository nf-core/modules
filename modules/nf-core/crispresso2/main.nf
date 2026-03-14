process CRISPRESSO2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crispresso2:2.3.3--py39hff726c5_0' :
        'biocontainers/crispresso2:2.3.3--py39hff726c5_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("CRISPResso_on_*")        , emit: results
    tuple val(meta), path("*.html")                 , emit: html
    tuple val(meta), path("CRISPResso_on_*/*.txt")  , emit: txt
    path "versions.yml"                             , emit: versions

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
    CRISPResso \\
        ${read_inputs} \\
        --name ${prefix} \\
        --output_folder . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crispresso2: \$(CRISPResso --version 2>&1 | grep -o 'CRISPResso version [0-9.]*' | sed 's/CRISPResso version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env bash
    set -e

    mkdir -p CRISPResso_on_${prefix}
    touch CRISPResso_on_${prefix}/${prefix}.html
    touch CRISPResso_on_${prefix}/CRISPResso_report.html
    touch CRISPResso_on_${prefix}/CRISPResso_mapping_statistics.txt
    touch CRISPResso_on_${prefix}/quantification_of_editing_frequency.txt
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crispresso2: 2.3.3
    END_VERSIONS
    """
}
