// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

params.options = [:]
options        = initOptions(params.options)

process CSVTK_SPLIT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(variants)
    val in_format
    val out_format

    output:
    tuple val(meta), path("*.${out_extension}"), emit: split_csv
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def delimiter = in_format == "tsv" ? "--tabs" : (in_format == "csv" ? "--delimiter ',' " : in_format)
    def out_delimiter = out_format == "tsv" ? "--out-tabs" : (out_format == "csv" ? "--out-delimiter ',' " : out_format)
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'

    """
    sed -i.bak '/^##/d' $variants
    csvtk \\
        split \\
        $options.args \\
        --num-cpus $task.cpus \\
        ${delimiter} \\
        ${out_delimiter} \\
        $variants


    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( csvtk version | sed -e 's/csvtk v//g' ))
    END_VERSIONS
    """
}
