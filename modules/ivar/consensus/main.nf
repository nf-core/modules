// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/ivar:1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*.fa")      , emit: fasta
    tuple val(meta), path("*.qual.txt"), emit: qual
    tuple val(meta), path("*.mpileup") , optional:true, emit: mpileup
    path "versions.yml"                , emit: versions

    script:
    def prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def save_mpileup = params.save_mpileup ? "tee ${prefix}.mpileup |" : ""
    """
    samtools mpileup \\
        --reference $fasta \\
        $options.args2 \\
        $bam | \\
        $save_mpileup \\
        ivar consensus \\
            $options.args \\
            -p $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """
}
