// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
process MINIA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::minia=3.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minia:3.2.4--he513fc3_0"
    } else {
        container "quay.io/biocontainers/minia:3.2.4--he513fc3_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fa'), emit: contigs
    tuple val(meta), path('*.unitigs.fa'), emit: unitigs
    tuple val(meta), path('*.h5')        , emit: h5
    path  '*.version.txt'                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "${reads.join("\n")}" > input_files.txt
    minia \\
        $options.args \\
        -nb-cores $task.cpus \\
        -in input_files.txt \\
        -out $prefix
    echo \$(minia --version 2>&1) | sed 's/^.*Minia version //; s/ .*\$//' > ${software}.version.txt
    """
}
