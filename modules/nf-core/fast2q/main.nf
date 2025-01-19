process FAST2Q {

    tag "2FAST2Q"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fast2q' :
        'biocontainers/2.7.2--pyh7e72e81_0' }"

    input:
    tuple val(meta1), path(fastq)
    tuple val(meta2), path(library)

    output:
    tuple val(meta), path("./2FAST2Q_output/*.*")   , emit: processed
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def input_file      = "--s ${fastq}" ?: ''
    def library_file    = library ? "--g ${library}" : ''

    """
    mkdir -p ./2FAST2Q_output

    2fast2q -c --o ./2FAST2Q_output/ $input_file $library_file $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
        2FAST2Q version: \$(2fast2q -v | grep 'Version:' | sed 's/Version: //g')
    END_VERSIONS
    """
}
