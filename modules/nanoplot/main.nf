// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NANOPLOT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::nanoplot=1.38.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanoplot:1.38.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png") , emit: png
    tuple val(meta), path("*.txt") , emit: txt
    tuple val(meta), path("*.log") , emit: log
    path  "versions.yml"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    def input_file = ("$ontfile".endsWith(".fastq.gz")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    """
    NanoPlot \\
        $options.args \\
        -t $task.cpus \\
        $input_file
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(NanoPlot --version 2>&1 | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
