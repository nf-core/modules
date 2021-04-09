// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SHOVILL {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::shovill=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/shovill:1.1.0--0"
    } else {
        container "quay.io/biocontainers/shovill:1.1.0--0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("contigs.fa")                         , emit: contigs
    tuple val(meta), path("shovill.corrections")                , emit: corrections
    tuple val(meta), path("shovill.log")                        , emit: log
    tuple val(meta), path("{skesa,spades,megahit,velvet}.fasta"), emit: raw_contigs
    tuple val(meta), path("contigs.{fastg,gfa,LastGraph}")      , optional:true, emit: gfa
    path "*.version.txt"                                        , emit: version

    script:
    def software = getSoftwareName(task.process)
    def memory = task.memory.toGiga()
    """
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $options.args \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force

    echo \$(shovill --version 2>&1) | sed 's/^.*shovill //'  > ${software}.version.txt
    """
}
