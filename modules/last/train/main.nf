// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1250--h2e03b76_0"
    }

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.par"), emit: param_file
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)

    last-train \\
        $options.args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        > ${prefix}.\$INDEX_NAME.par

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """
}
