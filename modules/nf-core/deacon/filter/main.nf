process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta), path(index), path(reads)

    output:
    tuple val(meta), path("${prefix}*.fq.gz"), emit: fastq_filtered
    tuple val(meta), path("${prefix}.json")  , emit: log
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_type = (reads instanceof List) ? "-o ${prefix}_1.fq -O ${prefix}_2.fq" : "> ${prefix}.fq" // deacon's automatic compression does not work
    """
    deacon \\
        filter \\
        --threads ${task.cpus} \\
        $args \\
        --summary ${prefix}.json \\
        -d $index \\
        $reads \\
        ${read_type}

    gzip -f ${prefix}*.fq
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > '${prefix}.fq.gz'
    touch ${prefix}.json
    """
}
