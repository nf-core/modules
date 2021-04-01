// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FLASH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    conda (params.enable_conda ? "bioconda::flash=1.2.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/flash:1.2.11--hed695b0_5"
    } else {
        container "quay.io/biocontainers/flash:1.2.11--hed695b0_5"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.*.fastq.gz"), emit: reads
    path "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def merged   = "-o ${prefix}.merged"
    def input_reads = "${reads[0]} ${reads[1]}"
    """
    flash \\
        $options.args \\
        $merged \\
        -z \\
        $input_reads
    echo \$(flash --version) > ${software}.version.txt
    """
}
