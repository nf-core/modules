process FLASH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flash:1.2.11--h5bf99c6_6' :
        'quay.io/biocontainers/flash:1.2.11--h5bf99c6_6' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.extendedFrags.fastq.gz"), emit: merged
    tuple val(meta), path("${prefix}.notCombined_*.fastq.gz"), emit: notcombined
    tuple val(meta), path("${prefix}.hist")                  , emit: histogram
    tuple val("${task.process}"), val('flash'), eval('flash --version 2>&1 | head -1 | sed \'s/^.*FLASH v//; s/ .*$//\''), topic: versions, emit: versions_flash

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("$reads" == "${prefix}.extendedFrags.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$reads" == "${prefix}.notCombined_1.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if ("$reads" == "${prefix}.notCombined_2.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    flash \\
        $args \\
        -o ${prefix} \\
        -z \\
        ${reads[0]} \\
        ${reads[1]}

    """


    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    echo "" | gzip > ${prefix}.extendedFrags.fastq.gz
    echo "" | gzip > ${prefix}.notCombined_1.fastq.gz
    echo "" | gzip > ${prefix}.notCombined_2.fastq.gz
    touch ${prefix}.hist

    """
}
