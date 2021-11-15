// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BUSTOOLS_SORT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bustools=0.41.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bustools:0.41.0--h60f4f9f_0"
    } else {
        container "quay.io/biocontainers/bustools:0.41.0--h60f4f9f_0"
    }

    input:
    tuple val(meta), path(bus)

    output:
    tuple val(meta), path("*.bus"), emit: bus
    path  "versions.yml"          , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bustools sort \\
        -o ${prefix}.bus \\
        -t $task.cpus \\
        -m ${task.memory.toGiga()}G \\
        $options.args \\
        $bus

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bustools version | sed -e "s/bustools, version //g")
    END_VERSIONS
    """
}
