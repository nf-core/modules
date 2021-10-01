
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PICARD_SORTSAM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::picard=2.25.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picard:2.25.7--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/picard:2.25.7--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)
    val sort_order

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard SortSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        SortSam \\
        -Xmx${avail_mem}g \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.bam \\
        --SORT_ORDER $sort_order

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(picard SortSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
