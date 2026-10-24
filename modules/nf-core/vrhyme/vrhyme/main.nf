// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VRHYME_VRHYME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9d994695f60fba485a0e2679c84965fa6fc2787e5b24df2c2d3e2969c7702ad5/data'
:         'community.wave.seqera.io/library/vrhyme_numpy_python_scikit-learn_setuptools:02d8cf67999994aa' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("vRhyme_best_bins_fasta/")                , emit: bins
    tuple val(meta), path("**/vRhyme_best_bins.*.membership.tsv")   , emit: membership
    tuple val(meta), path("**/vRhyme_best_bins.*.summary.tsv")      , emit: summary
    tuple val("${task.process}"), val('vrhyme'), eval("vRhyme --version 2>&1 | sed 's/^.*vRhyme v//; s/Using.*\$//'"), emit: versions_vrhyme, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_input = fasta.toString().replaceAll(/\.gz$/, '')
    def gunzip      = fasta.getExtension() == "gz" ? "gunzip -c ${fasta} > ${fasta_input}" : ""
    def cleanup     = fasta.getExtension() == "gz" ? "rm ${fasta_input}" : ""
    """
    ${gunzip}

    vRhyme \\
        -i ${fasta_input} \\
        -r ${reads} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}

    mv $prefix/vRhyme_best_bins_fasta/ vRhyme_best_bins_fasta

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/vRhyme_best_bins.19.membership.tsv
    touch ${prefix}/vRhyme_best_bins.19.summary.tsv
    mkdir -p vRhyme_best_bins_fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_1.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_10.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_11.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_12.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_13.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_14.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_2.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_3.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_4.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_5.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_6.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_7.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_8.fasta
    touch vRhyme_best_bins_fasta/vRhyme_bin_9.fasta
    """
}
