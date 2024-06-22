process PBCCS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbccs:6.4.0--h9ee0642_0' :
        'biocontainers/pbccs:6.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(pbi)
    val chunk_num
    val chunk_on

    output:
    tuple val(meta), path("*.chunk*.bam")     , emit: bam
    tuple val(meta), path("*.chunk*.bam.pbi") , emit: pbi
    tuple val(meta), path("*.report.txt" )    , emit: report_txt
    tuple val(meta), path("*.report.json" )   , emit: report_json
    tuple val(meta), path("*.metrics.json.gz"), emit: metrics
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ccs \\
        $bam \\
        ${prefix}.chunk${chunk_num}.bam \\
        --report-file ${prefix}.chunk${chunk_num}.report.txt \\
        --report-json ${prefix}.chunk${chunk_num}.report.json \\
        --metrics-json ${prefix}.chunk${chunk_num}.metrics.json.gz \\
        --chunk $chunk_num/$chunk_on \\
        -j $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbccs: \$(echo \$(ccs --version 2>&1) | grep 'ccs' | sed 's/^.*ccs //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch dummy.chunk1.bam
    touch dummy.chunk1.bam.pbi
    touch dummy.report.txt
    touch dummy.report.json
    echo "test" > dummy.metrics.json
    gzip dummy.metrics.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbccs: \$(echo \$(ccs --version 2>&1) | grep 'ccs' | sed 's/^.*ccs //; s/ .*\$//')
    END_VERSIONS
    """
}
