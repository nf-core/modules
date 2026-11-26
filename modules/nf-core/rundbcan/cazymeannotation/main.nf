// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RUNDBCAN_CAZYMEANNOTATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/31670e6b98b85351cf6263be7519e4e93e513d76a06e0c8f2e82a930e67dbd17/data'
        : 'community.wave.seqera.io/library/dbcan_pandas:617c1839ee862102'}"

    input:
    tuple val(meta), path(input_raw_data)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}_overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_diamond.out"), emit: dbcandiamond_results
    tuple val("${task.process}"), val('rundbcan'), eval("run_dbcan version | sed 's/dbCAN version: //g'"), emit: versions_rundbcan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_dbcan CAZyme_annotation \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir . \\
        --threads ${task.cpus} \\
        ${args}

    mv overview.tsv ${prefix}_overview.tsv
    mv dbCAN_hmm_results.tsv ${prefix}_dbCAN_hmm_results.tsv
    mv dbCANsub_hmm_results.tsv ${prefix}_dbCANsub_hmm_results.tsv
    mv diamond.out ${prefix}_diamond.out
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_overview.tsv
    touch ${prefix}_dbCAN_hmm_results.tsv
    touch ${prefix}_dbCANsub_hmm_results.tsv
    touch ${prefix}_diamond.out
    """
}
