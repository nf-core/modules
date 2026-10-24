// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RASUSA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7c52048426ed5860a231ffe6fb4c40acb23121969d87a231d5e1c7c3cc8661c/data'
:         'community.wave.seqera.io/library/rasusa:0.3.0--e22d5c77ebcbec79' }"

    input:
    tuple val(meta), path(reads), val(genome_size)
    val   depth_cutoff

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val("${task.process}"), val('rasusa'), eval('rasusa --version 2>&1 | sed -e "s/rasusa //g"'), emit: versions_rasusa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output   = meta.single_end ? "--output ${prefix}.fastq.gz" : "--output ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz"
    """
    rasusa \\
        $args \\
        --coverage $depth_cutoff \\
        --genome-size $genome_size \\
        --input $reads \\
        $output
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_1.fastq.gz
    if [ "${meta.single_end}" == "false" ]; then
        echo "" | gzip > ${prefix}_2.fastq.gz
    fi
    """
}
