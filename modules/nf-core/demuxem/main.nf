
process DEMUXEM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'docker'  : 'quay.io/cumulus/demuxem' }"
    
    input:
    tuple val(meta), path(input_raw_gene_bc_matrices_h5), path(input_hto_csv_file)
    val output_name
    val generate_gender_plot
    
    output:
    tuple val(meta), path("*_demux.zarr"), emit: zarr
    tuple val(meta), path("*.out.demuxEM.zarr"), emit: out_zarr
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def generateGenderPlot = generate_gender_plot ? "--generateGenderPlot ${generate_gender_plot}" : ""
    """

        demuxEM $input_raw_gene_bc_matrices_h5 \\
        $input_hto_csv_file $output_name \\
        $args \\
        $generateGenderPlot

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
           echo \$(demuxEM --version  2>&1)
        END_VERSIONS
        
    """

    stub:
    
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot 
    """
    touch ${prefix}.best

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echo \$(demuxEM --version  2>&1)
    END_VERSIONS
    """

}
