// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1=h089eab3_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/ivar:1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.fa")      , emit: fasta
    tuple val(meta), path("*.qual.txt"), emit: qual
    tuple val(meta), path("*.mpileup") , emit: mpileup
    path "*.version.txt"               , emit: version

    script:
    def software     = getSoftwareName(task.process)
    def prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def save_mpileup = options.save_mpileup ? "tee ${prefix}.mpileup |" : ""
    def min_bq = options.min_bq ? "-Q ${options.min_bq}" : "-Q 0"
    def max_depth = options.max_depth ? "-d ${options.max_depth}" : "-d 0"
    def count_orphans = options.count_orphans ? "-A" : ""
    def disable_baq = options.disable_baq ? "-B" : ""
    def output_all_bases = options.output_all_bases ? "-aa" : ""
    """
    samtools mpileup \\
        --fasta-ref $fasta \\
        $output_all_bases \\
        $count_orphans \\
        $disable_baq \\
        $max_depth \\
        $min_bq \\
        $options.args2 \\
        $bam | \\
        $save_mpileup \\
        ivar consensus \\
        $options.args \\
        -p $prefix

    ivar version | head -n1 2>&1 | sed 's/^.*iVar version //g' > ${software}.version.txt
    """
}
