// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ISOSEQ3_CLUSTER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::isoseq3=3.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/isoseq3:3.4.0--0"
    } else {
        container "quay.io/biocontainers/isoseq3:3.4.0--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.clustered.bam"),     path("*.clustered.singletons.bam"),       emit: bam
    tuple val(meta), path("*.clustered.bam.pbi"), path("*.clustered.singletons.bam.pbi"),   emit: bam_index
    path("*.{clustered.cluster,clustered.cluster_report.csv,clustered.transcriptset.xml}"), emit: reports
    path("*.{clustered.bam.pbi,clustered.hq.bam,clustered.hq.bam.pbi,clustered.lq.bam,clustered.lq.bam.pbi}"), emit: other_bams
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    cluster_out = bam.toString().replaceAll(/.bam$/, '.clustered.bam')
    """
    isoseq3 \\
        cluster \\
        $bam \\
        $cluster_out \\
        --singletons \\
        $options.args

    echo \$(isoseq3 --version 2>&1) | grep -e 'commit' > ${software}.version.txt
    """
}
