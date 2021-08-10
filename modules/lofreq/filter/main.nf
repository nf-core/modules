// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LOFREQ_FILTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4"
    } else {
        container "quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4"
    }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    lofreq \\
        filter \\
        $options.args \\
        -i $vcf \\
        -o ${prefix}.vcf.gz

    echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//' > ${software}.version.txt
    """
}
