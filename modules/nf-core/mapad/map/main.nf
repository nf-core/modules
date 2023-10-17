process MAPAD_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mapad=0.42.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapad:0.42.1--hc9368f3_2':
        'biocontainers/mapad:0.42.1--hc9368f3_2' }"

    input:
    tuple val(meta) , path(reads) // Supports only single-end or merged paired-end data
    tuple val(meta2), path(index)
    val mismatch_parameter
    val double_stranded_library
    val five_prime_overhang
    val three_prime_overhang
    val deam_rate_double_stranded
    val deam_rate_single_stranded
    val indel_rate

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def library_preparation = double_stranded_library ? 'double_stranded' : 'single_stranded'

    """
    INDEX=`find -L ./ -name "*.tbw" | sed 's/\\.tbw\$//'`

    mapad \\
        map \\
        ${args} \\
        --threads ${task.cpus} \\
        --reads ${reads} \\
        --reference \${INDEX} \\
        --output ${prefix}.bam \\
        -p ${mismatch_parameter} \\
        --library ${library_preparation} \\
        -f ${five_prime_overhang} \\
        -t ${three_prime_overhang} \\
        -d ${deam_rate_double_stranded} \\
        -s ${deam_rate_single_stranded} \\
        -i ${indel_rate}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapad: \$(echo \$(mapad --version) | sed 's/^mapAD //' ))
    END_VERSIONS
    """
}
