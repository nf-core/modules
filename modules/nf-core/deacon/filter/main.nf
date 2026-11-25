// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f1270fe75323ac6d213bf66d72b9f292c5ff31f7f5ee12300d2bed2b858f5ab/data'
        : 'community.wave.seqera.io/library/deacon:0.20.0--09d7e8e5c3f9dc88'}"

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
    if (!(reads instanceof List) && "${reads}" == "${prefix}.fq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
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
