// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MUSCLE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::muscle=3.8.1551" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/muscle:3.8.1551--h7d875b9_6"
    } else {
        container "quay.io/biocontainers/muscle:3.8.1551--h7d875b9_6"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.afa") , optional: true, emit: aligned_fasta
    tuple val(meta), path("*.phyi"), optional: true, emit: phyi
    tuple val(meta), path("*.phys"), optional: true, emit: phys
    tuple val(meta), path("*.clw") , optional: true, emit: clustalw
    tuple val(meta), path("*.html"), optional: true, emit: html
    tuple val(meta), path("*.msf") , optional: true, emit: msf
    tuple val(meta), path("*.tree"), optional: true, emit: tree
    path "*.log"                                   , emit: log
    path "versions.yml"                            , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fasta_out   = options.args.contains('-fasta') ? "-fastaout ${prefix}_muscle_msa.afa" : ''
    def clw_out     = options.args.contains('-clw') ? "-clwout ${prefix}_muscle_msa.clw" : ''
    def msf_out     = options.args.contains('-msf') ? "-msfout ${prefix}_muscle_msa.msf" : ''
    def phys_out    = options.args.contains('-phys') ? "-physout ${prefix}_muscle_msa.phys" : ''
    def phyi_out    = options.args.contains('-phyi') ? "-phyiout ${prefix}_muscle_msa.phyi" : ''
    def html_out    = options.args.contains('-html') ? "-htmlout ${prefix}_muscle_msa.html" : ''
    def tree_out    = options.args.contains('-maketree') ? "-out ${prefix}_muscle_msa.tree" : ''

    """
    muscle \\
        $options.args \\
        -in $fasta \\
        $fasta_out \\
        $clw_out \\
        $msf_out \\
        $phys_out \\
        $phyi_out \\
        $html_out \\
        $tree_out \\
        -loga muscle_msa.log
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(muscle -version |  sed 's/^MUSCLE v//; s/by.*\$//')
    END_VERSIONS
    """
}
