process KMERGENIE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f4197eec51307131e6cb0170a7969eda60995b23942d050f7495dc4a530b118/data':
        'community.wave.seqera.io/library/kmergenie:1.7051--675dfe5a4c7ea92b' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_report.html"), emit: html
    tuple val(meta), path("*.histo")      , emit: histo
    tuple val(meta), path("*.dat")        ,emit: dat
    tuple val(meta), path("*.pdf")        ,emit: pdf
    tuple val("${task.process}"), val('kmergenie'), eval('kmergenie --version |& sed "1!d ; s/KmerGenie //"'), emit: versions_kmergenie, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_list = reads.join("\n")

    """
    echo "$read_list" > ${prefix}_reads.txt

    kmergenie \\
        $args \\
        -o ${prefix} \\
        -t $task.cpus \\
        ${prefix}_reads.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}*.histo
    touch ${prefix}*.pdf
    """
}
