// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOLYSE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::nanolyse=1.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanolyse:1.2.0--py_0"
    } else {
        container "quay.io/biocontainers/nanolyse:1.2.0--py_0"
    }

    input:
    tuple val(meta), path(fastq)
    path nanolyse_fasta

    output:
    tuple val(meta), path("*.fastq.gz")   ,emit: nanolyse_fastq
    path "*.log"                          ,emit: nanolyse_log
    path "*.version.txt"                  ,emit: version

    script:
    """
    gunzip -c $fastq | NanoLyse -r $nanolyse_fasta | gzip > ${meta.id}.clean.fastq.gz
    cp NanoLyse.log ${meta.id}.nanolyse.log
    NanoLyse --version &> nanolyse.version.txt
    """
}
