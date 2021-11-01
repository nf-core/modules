// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUNC_DOWNLOADDB {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gunc=1.0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gunc:1.0.5--pyhdfd78af_0"
    }

    input:
    val db_name

    output:
    path "*.dmnd", emit: db
    path "versions.yml"          , emit: versions

    script:
    """
    gunc download_db . -db $db_name $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gunc --version )
    END_VERSIONS
    """
}
