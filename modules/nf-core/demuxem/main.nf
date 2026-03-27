process DEMUXEM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/demuxem:0.1.7.post1--pyhdfd78af_0' :
        'biocontainers/demuxem:0.1.7.post1--pyhdfd78af_0' }"
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
