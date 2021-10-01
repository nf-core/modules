// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0"
    } else {
        container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"
    }

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  index

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $reads \\
        -o ${prefix}.sam \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
