// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_MKREF {
    tag 'mkref'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ?: null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "nfcore/cellranger:6.0.2"
    } else {
        container "nfcore/cellranger:6.0.2"
    }

    input:
    path fasta
    path gtf
    val(reference_name)

    output:
    path "${reference_name}", emit: reference
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    cellranger mkref \\
        --genome=${reference_name} \\
        --fasta=${fasta} \\
        --genes=${gtf}

    cellranger --version | grep -o "[0-9\\. ]\\+" > ${software}.version.txt
    """
}
