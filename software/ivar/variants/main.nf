// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/ivar:1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.tsv"), emit: variants
    tuple val(meta), path("*.mpileup"), emit: mpileup
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def save_mpileup = params.containsKey('save_mpileup') ? "tee ${prefix}.mpileup |" : ""
    def gff      = params.containsKey('gff') ? "-g ${params.get('gff')}" : ""
    """
    samtools mpileup \\
        $options.args2 \\
        --reference $fasta \\
        $bam | \\
        $save_mpileup  \\
        ivar variants \\
        $options.args \\
        $gff \\
        -r $fasta \\
        -p $prefix

    ivar version | head -n1 2>&1 | sed 's/^.*iVar version //g' > ${software}.version.txt
    """
}
