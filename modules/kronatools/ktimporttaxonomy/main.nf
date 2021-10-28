// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process KRONATOOLS_KTIMPORTTAXONOMY {
    tag "${meta.classifier}-${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['classifier','id']) }

        conda (params.enable_conda ? "bioconda::krona=2.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/krona:2.8--pl5262hdfd78af_2"
    } else {
        container "quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2"
    }

    input:
    tuple val(meta), path(report)
    path  "taxonomy/taxonomy.tab"

    output:
    path("*.html")      , emit: html
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    ktImportTaxonomy "$report" -tax taxonomy
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(ktImportTaxonomy 2>&1) | sed 's/^.*KronaTools //; s/ - ktImportTaxonomy.*//')
    END_VERSIONS
    """
}
