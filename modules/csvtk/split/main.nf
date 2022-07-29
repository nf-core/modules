process CSVTK_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(csv)
    val in_format
    val out_format

    output:
    tuple val(meta), path("*.${out_extension}"), emit: split_csv
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def delimiter = in_format == "tsv" ? "--tabs" : (in_format == "csv" ? "--delimiter ',' " : in_format)
    def out_delimiter = out_format == "tsv" ? "--out-tabs" : (out_format == "csv" ? "--out-delimiter ',' " : out_format)
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    """
    sed -i.bak '/^##/d' $csv
    csvtk \\
        split \\
        $args \\
        --num-cpus $task.cpus \\
        $delimiter \\
        $out_delimiter \\
        $csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e 's/csvtk v//g' ))
    END_VERSIONS
    """
}
