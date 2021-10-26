// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CELLRANGER_MKREF {
    tag 'mkref'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
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
    path "versions.yml"     , emit: versions
    path "${reference_name}", emit: reference

    script:
    """
    cellranger mkref \\
        --genome=${reference_name} \\
        --fasta=${fasta} \\
        --genes=${gtf}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | grep -o "[0-9\\. ]\\+" )
    END_VERSIONS
    """
}
