
process DEMUXEM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
        container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'cumulusprod/demuxem:0.1.7' }"
    input:

    tuple val(meta), path(input_raw_gene_bc_matrices_h5), path(input_hto_csv_file)
    val output_name
    val generate_gender_plot
    when:
    task.ext.when == null || task.ext.when

    output:
    tuple val(meta), path("*_demux.zarr"),path("*.out.demuxEM.zarr"), emit: results
    path "versions.yml"                        , emit: versions
    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def generateGenderPlot = generate_gender_plot != 'None' ? " --generateGenderPlot 
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
