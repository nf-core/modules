// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEMUXEM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7c6c63164c027a83143685618eba0cae8b814c6ffe87c2379b7c9b059402fd0/data'
        : 'community.wave.seqera.io/library/demuxem_numpy_pandas_python:24f4768270179af6'}"
    input:
    tuple val(meta), path(input_raw_gene_bc_matrices_h5), path(input_hto_csv_file)
    val generate_gender_plot
    val genome
    val generate_diagnostic_plots

    output:
    tuple val(meta), path("*_demux.zarr.zip"), emit: zarr
    tuple val(meta), path("*.out.demuxEM.zarr.zip"), emit: out_zarr
    tuple val("${task.process}"), val('demuxEM'), eval("demuxEM --version"), topic: versions, emit: versions_demuxem

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def generateGenderPlot = generate_gender_plot ? "--generate-gender-plot $generate_gender_plot" : ""
    def genome_file        = genome ? "--genome $genome" : ""
    def diagnostic_plots   = generate_diagnostic_plots ? "--generate-diagnostic-plots $generate_diagnostic_plots" : ""
    """
    demuxEM \\
        ${args} \\
        ${generateGenderPlot}\\
        ${genome_file}\\
        ${diagnostic_plots} \\
        ${input_raw_gene_bc_matrices_h5} \\
        ${input_hto_csv_file} \\
        ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out.demuxEM.zarr.zip
    touch ${prefix}_demux.zarr.zip
    """

}
