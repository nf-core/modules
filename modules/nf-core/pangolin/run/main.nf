// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PANGOLIN_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e718928024db0f5501a189533e70779345e8d3e0c7782b09e9c459d665ab15a/data' :
        'community.wave.seqera.io/library/pangolin-data_pangolin_snakemake-minimal:a1e2ba4975d2f9f2' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.csv'), emit: report
    tuple val("${task.process}"), val('pangolin'), eval("pangolin -v | sed 's/pangolin //'"), emit: versions_pangolin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_command = db ? "--datadir ${db}" : ''
    """
    export XDG_CACHE_HOME=/tmp/.cache

    pangolin \\
        $fasta\\
        $db_command \\
        --outfile ${prefix}.pangolin.csv \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=/tmp/.cache

    touch ${prefix}.pangolin.csv
    """
}
