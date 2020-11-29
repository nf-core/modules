// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options = initOptions(params.options)

process MULTIQC {
    tag "multiqc"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id:'') }

    conda(params.enable_conda ? "bioconda::multiqc=1.9" : null)
    container "quay.io/biocontainers/multiqc:1.9--py_1"

    input:
    path(generic_report)

    output:
    path("multiqc_data"), emit: dir
    path("multiqc_report.html"), emit: html
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    multiqc .
    multiqc --version | sed -e "s/multiqc, version //g" > ${software}.version.txt
    """
}
