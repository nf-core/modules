// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ISOSEQ3_REFINE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::isoseq3=3.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/isoseq3:3.4.0--0"
    } else {
        container "quay.io/biocontainers/isoseq3:3.4.0--0"
    }

    input:
    tuple val(meta), path(bam), path(primers)

    output:
    tuple val(meta), path("*.flnc.bam"), path("*.flnc.bam.pbi"), emit: bam
    tuple path("*{.consensusreadset.xml,filter_summary.json,.report.csv}"), emit: reports
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    refine_out = bam.toString().replaceAll(/.bam$/, '.flnc.bam')
    """
    isoseq3 \\
        refine \\
        -j $task.cpus \\
        $options.args \\
        $bam \\
        $primers \\
        $refine_out

    echo \$(isoseq3 refine --version) | grep -e 'commit' > ${software}.version.txt
    """
}
