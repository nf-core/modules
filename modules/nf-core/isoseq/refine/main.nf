process ISOSEQ_REFINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:4.0.0--h9ee0642_0' :
        'quay.io/biocontainers/isoseq:4.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)
    path primers

    output:
    tuple val(meta), path("*.bam")                       , emit: bam
    tuple val(meta), path("*.bam.pbi")                   , emit: pbi
    tuple val(meta), path("*.consensusreadset.xml")      , emit: consensusreadset
    tuple val(meta), path("*.filter_summary.report.json"), emit: summary
    tuple val(meta), path("*.report.csv")                , emit: report
    tuple val("${task.process}"), val('isoseq'), eval("isoseq refine --version | head -n 1 | sed 's/isoseq refine //' | sed 's/ (commit.\\+//'"), emit: versions_isoseq, topic: versions


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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.pbi
    touch ${prefix}.consensusreadset.xml
    touch ${prefix}.filter_summary.report.json
    touch ${prefix}.report.csv
    """
}
