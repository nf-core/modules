// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CLONALFRAMEML {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::clonalframeml=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/clonalframeml:1.12--h7d875b9_1"
    } else {
        container "quay.io/biocontainers/clonalframeml:1.12--h7d875b9_1"
    }

    input:
    tuple val(meta), path(newick), path(msa)

    output:
    tuple val(meta), path("*.emsim.txt")                   , emit: emsim, optional: true
    tuple val(meta), path("*.em.txt")                      , emit: em
    tuple val(meta), path("*.importation_status.txt")      , emit: status
    tuple val(meta), path("*.labelled_tree.newick")        , emit: newick
    tuple val(meta), path("*.ML_sequence.fasta")           , emit: fasta
    tuple val(meta), path("*.position_cross_reference.txt"), emit: pos_ref
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ClonalFrameML \\
        $newick \\
        <(gzip -cdf $msa) \\
        $prefix \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(ClonalFrameML -version 2>&1) | sed 's/^.*ClonalFrameML v//' )
    END_VERSIONS
    """
}
