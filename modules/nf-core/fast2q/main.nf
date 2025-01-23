process FAST2Q {

    tag "2FAST2Q"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fast2q' :
        'biocontainers/fast2q:2.7.2--pyh7e72e81_0' }"

    input:
    tuple val(meta1), path(fastq)
    tuple val(meta2), path(library)

    output:
    tuple val(meta1), path("2FAST2Q_output_*")   , emit: out_folder
    path "2FAST2Q_output_*/compiled.csv"         , emit: count_matrix
    path "2FAST2Q_output_*/compiled_stats.csv"   , emit: run_statistics
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def input_file      = "--s ${fastq}" ?: ''
    def library_file    = (library instanceof Path && library.exists()) ? "--g ${library}" : ''

    """

    2fast2q -c --o ./ $input_file $library_file $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
        2FAST2Q version: \$(2fast2q -v | grep 'Version:' | sed 's/Version: //g')
    END_VERSIONS
    """
}