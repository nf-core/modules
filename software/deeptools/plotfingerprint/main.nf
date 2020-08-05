// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process DEEPTOOLS_PLOTFINGERPRINT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/deeptools:3.4.3--py_0"
    //container "https://depot.galaxyproject.org/singularity/deeptools:3.4.3--py_0"

    conda (params.conda ? "bioconda::deeptools=3.4.3" : null)

    input:
    tuple val(meta), path(bams), path(bais)
    val options

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.raw.txt"), emit: matrix
    tuple val(meta), path("*.qcmetrics.txt"), emit: metrics
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def extend   = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    """
    plotFingerprint \\
        $ioptions.args \\
        $extend \\
        --bamfiles ${bams.join(' ')} \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors $task.cpus

    plotFingerprint --version | sed -e "s/plotFingerprint //g" > ${software}.version.txt
    """
}
