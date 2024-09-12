process ISOSEQ_REFINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0' :
        'biocontainers/isoseq:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    path primers

    output:
    tuple val(meta), path("*.bam")                       , emit: bam
    tuple val(meta), path("*.bam.pbi")                   , emit: pbi
    tuple val(meta), path("*.consensusreadset.xml")      , emit: consensusreadset
    tuple val(meta), path("*.filter_summary.report.json"), emit: summary
    tuple val(meta), path("*.report.csv")                , emit: report
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isoseq \\
        refine \\
        -j $task.cpus \\
        $args \\
        $bam \\
        $primers \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq refine --version | head -n 1 | sed 's/isoseq refine //' | sed 's/ (commit.\\+//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch dummy.bam
    touch dummy.bam.pbi
    touch dummy.consensusreadset.xml
    touch dummy.filter_summary.report.json
    touch dummy.report.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoseq: \$( isoseq refine --version | head -n 1 | sed 's/isoseq refine //' | sed 's/ (commit.\\+//' )
    END_VERSIONS
    """
}
