process UCSC_GTFTOGENEPRED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:482--h0b57e2e_0':
        'biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.genepred"), emit: genepred
    tuple val(meta), path("*.refflat") , emit: refflat , optional: true
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gen_refflat = args.contains("-genePredExt") && args.contains("-geneNameAsName2") ? "true" : "false"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gtfToGenePred \\
        $args \\
        $gtf  \\
        ${prefix}.genepred

    if [ "${gen_refflat}" == "true" ] ; then
        awk 'BEGIN { OFS="\\t"} {print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' ${prefix}.genepred > ${prefix}.refflat
    fi

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genepred
    touch ${prefix}.refflat

    """
}
