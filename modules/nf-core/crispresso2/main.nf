process CRISPRESSO2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crispresso2:2.3.3--py39hff726c5_0' :
        'biocontainers/crispresso2:2.3.3--py39hff726c5_0' }"

    input:
    tuple val(meta), path(reads)
    val amplicon_seq_fallback
    val guide_seq_fallback

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
    
    // Use amplicon and guide from meta if available, otherwise fallback to input parameters
    def amplicon_seq = meta.amplicon_seq ?: amplicon_seq_fallback ?: ""
    def guide_seq = meta.guide_seq ?: guide_seq_fallback ?: ""
    
    def amplicon_param = amplicon_seq ? "-a ${amplicon_seq}" : ""
    def guide_param = guide_seq ? "-g ${guide_seq}" : ""
    
    // Handle single-end vs paired-end reads and copy files to avoid symlink issues in Docker
    def read_inputs = ""
    def copy_commands = ""
    if (reads instanceof Path || reads.size() == 1) {
        copy_commands = "cp -L ${reads} reads.fastq.gz"
        read_inputs = "-r1 reads.fastq.gz"
    } else {
        copy_commands = "cp -L ${reads[0]} reads_1.fastq.gz\n    cp -L ${reads[1]} reads_2.fastq.gz"
        read_inputs = "-r1 reads_1.fastq.gz -r2 reads_2.fastq.gz"
    }

    """
    ${copy_commands}
    
    CRISPResso \\
        ${read_inputs} \\
        ${amplicon_param} \\
        ${guide_param} \\
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
    // Ensure stub runs without Docker interference
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
