// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

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
    tuple val(meta), path("*.bam")               , emit: bam
    tuple val(meta), path("*.bam.pbi")           , emit: pbi
    tuple val(meta), path("*.cluster")           , emit: cluster
    tuple val(meta), path("*.cluster_report.csv"), emit: cluster_report
    tuple val(meta), path("*.transcriptset.xml") , emit: transcriptset
    tuple val(meta), path("*.hq.bam")            , emit: hq_bam
    tuple val(meta), path("*.hq.bam.pbi")        , emit: hq_pbi
    tuple val(meta), path("*.lq.bam")            , emit: lq_bam
    tuple val(meta), path("*.lq.bam.pbi")        , emit: lq_pbi
    path  "versions.yml"                         , emit: versions

    tuple val(meta), path("*.singletons.bam")    , optional: true, emit: singletons_bam
    tuple val(meta), path("*.singletons.bam.pbi"), optional: true, emit: singletons_pbi

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    isoseq3 \\
        cluster \\
        $bam \\
        ${prefix}.bam \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        isoseq3 cluster: \$( isoseq3 cluster --version|sed 's/isoseq cluster //g'|sed 's/ (.*//g' )
    END_VERSIONS
    """
}
