process PCNE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pcne%3A3.3.1--hdfd78af_0'
        : 'quay.io/biocontainers/pcne:3.3.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(chromosome)
    tuple val(meta3), path(plasmids)

    output:
    tuple val(meta), path("*_results.tsv"), emit: results
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.png"), emit: plots, optional: true
    tuple val("${task.process}"), val('pcne'), eval("pcne -v | grep 'version' | tail -n 1 | sed 's/.*version //'"), topic: versions, emit: versions_pcne

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def is_bam = reads.toString().endsWith('.bam')
    def input_cmd = is_bam ? "-b ${reads}" : "-r ${reads[0]} -R ${reads[1]}"

    def plasmid_cmd = plasmids instanceof List ? plasmids.join(' ') : plasmids

    """
    pcne \\
        -c ${chromosome} \\
        -p ${plasmid_cmd} \\
        ${input_cmd} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.tsv
    touch ${prefix}.log
    """
}
