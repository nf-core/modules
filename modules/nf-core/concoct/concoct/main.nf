// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CONCOCT_CONCOCT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b414706a0ce09cd26055f920ea4d7db8e6e0a4845f806f0918e9792a2c847f95/data'
        : 'community.wave.seqera.io/library/concoct_setuptools:7d62fa8e6704b1fe'}"

    input:
    tuple val(meta), path(coverage_file), path(fasta)

    output:
    tuple val(meta), path("*_args.txt"), emit: args_txt
    tuple val(meta), path("*_clustering_gt1000.csv"), emit: clustering_csv
    tuple val(meta), path("*_log.txt"), emit: log_txt
    tuple val(meta), path("*_original_data_gt1000.csv"), emit: original_data_csv
    tuple val(meta), path("*_PCA_components_data_gt1000.csv"), emit: pca_components_csv
    tuple val(meta), path("*_PCA_transformed_data_gt1000.csv"), emit: pca_transformed_csv
    tuple val("${task.process}"), val('concoct'), eval('concoct --version | cut -d " " -f2'), emit: versions_concoct, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    concoct \\
        ${args} \\
        --threads ${task.cpus} \\
        --coverage_file ${coverage_file} \\
        --composition_file ${fasta} \\
        -b ${prefix}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_args.txt
    touch ${prefix}_clustering_gt1000.csv
    touch ${prefix}_log.txt
    touch ${prefix}_original_data_gt1000.csv
    touch ${prefix}_PCA_components_data_gt1000.csv
    touch ${prefix}_PCA_transformed_data_gt1000.csv

    """
}
