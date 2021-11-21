// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PEDDY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::peddy=0.4.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/peddy:0.4.8--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/peddy:0.4.8--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path ped

    output:
    tuple val(meta), path("*.html")     , emit: html
    tuple val(meta), path("*.csv")      , emit: csv
    tuple val(meta), path("*.peddy.ped"), emit: ped
    tuple val(meta), path("*.png")      , emit: png
    path "versions.yml"                 , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    peddy \\
        $options.args \\
        --plot \\
        -p $task.cpus \\
        $vcf \\
        $ped

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( peddy --version 2>&1 | sed 's/peddy, version //' )
    END_VERSIONS
    """
}
