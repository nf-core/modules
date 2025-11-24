process CSVTK_CONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.31.0--h9ee0642_0' :
        'biocontainers/csvtk:0.31.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(csv, name: 'inputs/csv*/*')
    val in_format
    val out_format

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: csv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: (out_format == 'tsv' ? 'tsv' : 'csv')
    def delimiter = in_format == "tsv" ? "\t" : (in_format == "csv" ? "," : in_format)
    def out_delimiter = out_format == "tsv" ? "\t" : (out_format == "csv" ? "," : out_format)
    """
    csvtk \\
        concat \\
        $args \\
        --num-cpus $task.cpus \\
        --delimiter "${delimiter}" \\
        --out-delimiter "${out_delimiter}" \\
        --out-file ${prefix}.${suffix} \\
        $csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: (out_format == 'tsv' ? 'tsv' : 'csv')
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
