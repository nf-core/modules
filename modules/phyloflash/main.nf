// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PHYLOFLASH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::phyloflash=3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/phyloflash:3.4--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/phyloflash:3.4--hdfd78af_1"
    }

    input:
    tuple val(meta), path(reads)
    path(database)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    //TODO Accomodate interleaved
    """
    phyloFlash.pl \\
        ${options.args} \\
        -read1 ${reads[0]} \\
        -read2 ${reads[1]} \\
        -CPUs ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "sed 's/^.*phyloFlash //")
    END_VERSIONS
    """
}
