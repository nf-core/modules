// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3"
    } else {
        container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    }

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "$reads" : "${reads[0]} ${reads[1]}"
    """
    minimap2 \\
        $options.args \\
        -t $task.cpus \\
        $reference \\
        $input_reads \\
        > ${prefix}.paf

    echo \$(minimap2 --version 2>&1) > ${software}.version.txt
    """
}
