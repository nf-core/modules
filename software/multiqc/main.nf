// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options = initOptions(params.options)

process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: meta.id) }

    conda(params.enable_conda ? "bioconda::multiqc=1.9" : null)
    container "quay.io/biocontainers/multiqc:1.9--py_1"

    input:
    tuple val(meta), path("*.html")

    output:
    tuple val(meta), path("multiqc_data"), emit: dir
    tuple val(meta), path("*.html"), emit: html
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
        """
        multiqc --version | sed -e "s/version//g" > ${software}.version.txt
        """
    }
}
