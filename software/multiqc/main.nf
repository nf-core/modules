// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    tag "multiqc"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    conda(params.enable_conda ? "bioconda::multiqc=1.9" : null)
    container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"

    input:
//    path multiqc_config
//    path multiqc_custom_config
    path('fastqc/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "*_plots", optional: true, emit: plots

    script:
//    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $options.args $custom_config .
    """
}
