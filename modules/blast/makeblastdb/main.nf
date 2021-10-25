// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0'
    } else {
        container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'
    }

    input:
    path fasta

    output:
    path 'blast_db'     , emit: db
    path "versions.yml" , emit: versions

    script:
    """
    makeblastdb \\
        -in $fasta \\
        $options.args
    mkdir blast_db
    mv ${fasta}* blast_db
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
