process ISOSEQ_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0' :
        'quay.io/biocontainers/isoseq:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.transcripts.bam")               , emit: bam
    tuple val(meta), path("*.transcripts.bam.pbi")           , emit: pbi
    tuple val(meta), path("*.transcripts.cluster")           , emit: cluster
    tuple val(meta), path("*.transcripts.cluster_report.csv"), emit: cluster_report
    tuple val(meta), path("*.transcripts.transcriptset.xml") , emit: transcriptset
    tuple val(meta), path("*.transcripts.hq.bam")            , optional: true, emit: hq_bam
    tuple val(meta), path("*.transcripts.hq.bam.pbi")        , optional: true, emit: hq_pbi
    tuple val(meta), path("*.transcripts.lq.bam")            , optional: true, emit: lq_bam
    tuple val(meta), path("*.transcripts.lq.bam.pbi")        , optional: true, emit: lq_pbi
    tuple val(meta), path("*.transcripts.singletons.bam")    , optional: true, emit: singletons_bam
    tuple val(meta), path("*.transcripts.singletons.bam.pbi"), optional: true, emit: singletons_pbi
    tuple val("${task.process}"), val('isoseq'), eval("isoseq cluster --version | head -n 1 | sed 's/isoseq cluster //g' | sed 's/ (.*//g'"), emit: versions_isoseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/isoseq/cluster2

    Reason:
    isoseq cluster is being phased out upstream in favour of cluster2, which
    is faster and has no memory constraints on large read counts.
    """
    assert false: deprecation_message

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/isoseq/cluster2

    Reason:
    isoseq cluster is being phased out upstream in favour of cluster2, which
    is faster and has no memory constraints on large read counts.
    """
    assert false: deprecation_message
}
