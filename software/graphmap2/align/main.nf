// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda     (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"

    input:
    tuple val(meta), path(fastq)
    path(fasta)
    path(index)

    output:
    tuple val(meta), path("*.sam"), emit: align_sam

    script:
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $fastq \\
        -o ${meta.id}.sam \\
        --extcigar
    """
}
