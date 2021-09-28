// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NEXTCLADE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::nextclade_js=0.14.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nextclade_js:0.14.4--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/nextclade_js:0.14.4--h9ee0642_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.csv")       , emit: csv
    tuple val(meta), path("${prefix}.json")      , emit: json
    tuple val(meta), path("${prefix}.tree.json") , emit: json_tree
    tuple val(meta), path("${prefix}.tsv")       , emit: tsv
    tuple val(meta), path("${prefix}.clades.tsv"), optional:true, emit: tsv_clades
    path "versions.yml"                          , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    nextclade \\
        $options.args \\
        --jobs $task.cpus \\
        --input-fasta $fasta \\
        --output-json ${prefix}.json \\
        --output-csv ${prefix}.csv \\
        --output-tsv ${prefix}.tsv \\
        --output-tsv-clades-only ${prefix}.clades.tsv \\
        --output-tree ${prefix}.tree.json

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(nextclade --version 2>&1)
    END_VERSIONS
    """
}
