// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRINGTIE_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.4--h7e0af3c_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"
    }

    input:
    path stringtie_gtf
    path annotation_gtf

    output:
    path "stringtie.merged.gtf", emit: gtf

    script:
    """
    stringtie \\
        --merge $stringtie_gtf \\
        -G $annotation_gtf \\
        -o stringtie.merged.gtf
    """
}
