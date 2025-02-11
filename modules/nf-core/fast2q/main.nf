process FAST2Q {

    tag "2FAST2Q"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?

        'https://depot.galaxyproject.org/singularity/fast2q:2.7.4--pyh7e72e81_0' :
        'biocontainers/fast2q:2.7.4--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(library)

    output:
    tuple val(meta), path("${prefix}.csv")                      , emit: count_matrix
    tuple val(meta), path("${prefix}_stats.csv")                , emit: stats
    tuple val(meta), path("${prefix}_distribution_plot.png")    , emit: distribution_plot
    tuple val(meta), path("${prefix}_reads_plot.png")           , emit: reads_plot
    tuple val(meta), path("${prefix}_reads_plot_percentage.png"), emit: reads_plot_percentage
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    prefix              = task.ext.prefix ?: "${meta.id}"
    def input_file      = (fastq instanceof Path && fastq.exists()) ? "--s ${fastq}" : ''
    def library_file    = (library instanceof Path && library.exists()) ? "--g ${library}" : ''

    """
    export MPLCONFIGDIR=\$PWD
    2fast2q \\
        -c \\
        --o ./ \\
        --fn ${prefix} \\
        --cp ${task.cpus} \\
        $input_file \\
        $library_file \\
        $args

    mv **/${prefix}* .

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        2FAST2Q version: \$(2fast2q -v | grep 'Version:' | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    prefix              = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.csv
    touch ${prefix}_stats.csv
    touch ${prefix}_distribution_plot.png
    touch ${prefix}_reads_plot.png
    touch ${prefix}_reads_plot_percentage.png

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        2FAST2Q version: \$(2fast2q -v | grep 'Version:' | sed 's/Version: //g')
    END_VERSIONS
    """

}
