process ISOSEQ3_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::isoseq3=3.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/isoseq3:3.4.0--0"
    } else {
        container "quay.io/biocontainers/isoseq3:3.4.0--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.transcripts.bam")               , emit: bam
    tuple val(meta), path("*.transcripts.bam.pbi")           , emit: pbi
    tuple val(meta), path("*.transcripts.cluster")           , emit: cluster
    tuple val(meta), path("*.transcripts.cluster_report.csv"), emit: cluster_report
    tuple val(meta), path("*.transcripts.transcriptset.xml") , emit: transcriptset
    path  "versions.yml"                                     , emit: versions

    tuple val(meta), path("*.transcripts.hq.bam")            , optional: true, emit: hq_bam
    tuple val(meta), path("*.transcripts.hq.bam.pbi")        , optional: true, emit: hq_pbi
    tuple val(meta), path("*.transcripts.lq.bam")            , optional: true, emit: lq_bam
    tuple val(meta), path("*.transcripts.lq.bam.pbi")        , optional: true, emit: lq_pbi
    tuple val(meta), path("*.transcripts.singletons.bam")    , optional: true, emit: singletons_bam
    tuple val(meta), path("*.transcripts.singletons.bam.pbi"), optional: true, emit: singletons_pbi


    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    isoseq3 \\
        cluster \\
        $bam \\
        ${prefix}.transcripts.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        isoseq3 cluster: \$( isoseq3 cluster --version|sed 's/isoseq cluster //g'|sed 's/ (.*//g' )
    END_VERSIONS
    """
}
