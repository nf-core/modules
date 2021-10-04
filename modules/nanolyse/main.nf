// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

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
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.log"                       , emit: log
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gunzip -c $fastq | NanoLyse -r $fasta | gzip > ${prefix}.fastq.gz
    mv NanoLyse.log ${prefix}.nanolyse.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(NanoLyse --version 2>&1 | sed -e "s/NanoLyse //g")
    END_VERSIONS
    """
}
