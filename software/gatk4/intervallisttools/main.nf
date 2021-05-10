// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_INTERVALLISTTOOLS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(interval_list)

    output:
    tuple val(meta), path("*.split/*/*.interval_list"), emit: interval_list
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def out_dir = ${prefix}.split
    """

    mkdir ${out_dir}

    gatk \\
    --java-options "-Xms1g" \\
    IntervalListTools \\
    -I ${interval_list} \\
    -O ${out_dir} \\
    $options.args

    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
