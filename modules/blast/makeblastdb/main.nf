// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    path fasta

    output:
    path 'blast_db'     , emit: db
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    makeblastdb \\
        -in $fasta \\
        $options.args
    mkdir blast_db
    mv ${fasta}* blast_db
    echo \$(blastn -version 2>&1) | sed 's/^.*blastn: //; s/ .*\$//' > ${software}.version.txt
    """
}
