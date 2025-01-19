process FAST2Q {

    tag "2FAST2Q"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fast2q' :
        'quay.io/biocontainers/fast2q' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("./2FAST2Q_output/*.*")   , emit: processed
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def input   = "--s ${fastq}" ?: ''

    """
    mkdir ./2FAST2Q_output

    2fast2q -c --o ./2FAST2Q_output/ $input $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | sed 's/Python //g')
        2FAST2Q version: \$(2fast2q -v | grep 'Version:' | sed 's/Version: //g')
    END_VERSIONS
    """
}
