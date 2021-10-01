// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    tuple val(meta), path(bam)
    path index

    output:
    tuple val(meta), path("*.bedGraph.gz")         , emit: bedgraph
    tuple val(meta), path("*.txt.gz")              , emit: methylation_calls
    tuple val(meta), path("*.cov.gz")              , emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt")          , emit: mbias
    path "versions.yml"                            , emit: version

    script:
    def seqtype  = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --report \\
        $seqtype \\
        $options.args \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
