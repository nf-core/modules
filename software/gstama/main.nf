// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        //TODO Update singularity image path
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)
    path genome

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    tama_collapse.py \\
        -s $bam \\
        -f $genome \\
        -p $prefix \\
        $options.args \\

    echo \$(tama_collapse.py -version 2>&1) | grep 'tc_version_date_' > ${software}.version.txt
    """
}
