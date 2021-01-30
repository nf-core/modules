// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::methyldackel=0.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/methyldackel:0.5.0--hed50d52_0"
    } else {
        container "quay.io/biocontainers/methyldackel:0.5.0--hed50d52_0"
    }

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.txt") , emit: mbias
    path  "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    MethylDackel mbias \\
        $options.args \\
        $fasta \\
        $bam \\
        $prefix \\
        --txt \\
        > ${prefix}.txt

    echo \$(methyldackel --version 2>&1) | cut -f1 -d" " > ${software}.version.txt
    """
}
