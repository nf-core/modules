// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KRONATOOLS_KRONADB {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::krona=2.7.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/krona:2.7.1--pl526_5"
    } else {
        container "quay.io/biocontainers/krona:2.7.1--pl526_5"
    }
    input:

    output:
    path 'taxonomy/taxonomy.tab', emit: db
    path "versions.yml"         , emit: versions

    script:
    def VERSION='2.7.1'
    """
    ktUpdateTaxonomy.sh taxonomy

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: $VERSION
    END_VERSIONS
    """
}
