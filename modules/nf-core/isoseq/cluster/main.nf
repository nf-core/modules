process ISOSEQ_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0' :
        'biocontainers/isoseq:4.0.0--h9ee0642_0' }"

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
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isoseq \\
        cluster \\
        $bam \\
        ${prefix}.transcripts.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq cluster --version | head -n 1 | sed 's/isoseq cluster //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch dummy.transcripts.bam
    touch dummy.transcripts.bam.pbi
    touch dummy.transcripts.cluster
    touch dummy.transcripts.cluster_report.csv
    touch dummy.transcripts.transcriptset.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq cluster --version | head -n 1 | sed 's/isoseq cluster //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
