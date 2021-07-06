// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_LASTAL {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::last=1238" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1238--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1238--h2e03b76_0"
    }

    input:
    tuple val(meta), path(fastx), path (param_file)
    path index

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def trained_params = param_file ? "-p ${param_file}"  : ''
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)
    lastal \\
        $trained_params \\
        $options.args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        | gzip --no-name > ${prefix}.\$INDEX_NAME.maf.gz
    # gzip needs --no-name otherwise it puts a timestamp in the file,
    # which makes its checksum non-reproducible.

    echo \$(lastal --version 2>&1) | sed 's/lastal //' > ${software}.version.txt
    """
}
