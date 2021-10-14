// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEEPTOOLS_PLOTFINGERPRINT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0"
    } else {
        container "quay.io/biocontainers/deeptools:3.5.1--py_0"
    }

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("*.pdf")          , emit: pdf
    tuple val(meta), path("*.raw.txt")      , emit: matrix
    tuple val(meta), path("*.qcmetrics.txt"), emit: metrics
    path  "versions.yml"                    , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def extend   = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : ''
    """
    plotFingerprint \\
        $options.args \\
        $extend \\
        --bamfiles ${bams.join(' ')} \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors $task.cpus

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(plotFingerprint --version | sed -e "s/plotFingerprint //g")
    END_VERSIONS
    """
}
