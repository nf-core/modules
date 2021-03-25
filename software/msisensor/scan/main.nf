// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSOR_SCAN {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2"
    } else {
        container "quay.io/biocontainers/msisensor:0.5--hb3646a4_2"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple (val(meta), path("*.tab"), emit: txt)
    path ("*.version.txt"          , emit: version)

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    msisensor \\
        scan \\
        -d $fasta \\
        -o ${prefix}.msisensor_scan.tab \\
        $options.args

    echo \$(msisensor 2>&1) | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p' > ${software}.version.txt
    """
}
