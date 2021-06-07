// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NEXTCLADE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::nextclade_js=0.14.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nextclade_js:0.14.2--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/nextclade_js:0.14.2--h9ee0642_0"
    }

    input:
    tuple val(meta), path(fasta)
    val   output_format

    output:
    tuple val(meta), path("${prefix}.csv")       , optional:true, emit: csv
    tuple val(meta), path("${prefix}.json")      , optional:true, emit: json
    tuple val(meta), path("${prefix}.tree.json") , optional:true, emit: json_tree
    tuple val(meta), path("${prefix}.tsv")       , optional:true, emit: tsv
    tuple val(meta), path("${prefix}.clades.tsv"), optional:true, emit: tsv_clades
    path "*.version.txt"                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def format   = output_format
    if (!(format in ['json', 'csv', 'tsv', 'tree', 'tsv-clades-only'])) {
        format = 'json'
    }
    def extension = format
    if (format in ['tsv-clades-only']) {
        extension = '.clades.tsv'
    } else if (format in ['tree']) {
        extension = 'tree.json'
    }
    """
    nextclade \\
        $options.args \\
        --jobs $task.cpus \\
        --input-fasta $fasta \\
        --output-${format} ${prefix}.${extension}

    echo \$(nextclade --version 2>&1) > ${software}.version.txt
    """
}
