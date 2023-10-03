process UCSC_GTFTOGENEPRED {
    tag '${meta.id}'
    label 'process_low'

    conda "bioconda::ucsc-gtftogenepred=447"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:447--h954228d_0':
        'biocontainers/ucsc-gtftogenepred:447--h954228d_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.genepred"), emit: genepred
    tuple val(meta), path("*.refflat") , emit: refflat , optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gen_refflat = args.contains("-genePredExt") && args.contains("-geneNameAsName2") ? "true" : "false"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gtfToGenePred \\
        $args \\
        $gtf  \\
        ${prefix}.genepred

    if [ "${gen_refflat}" == "true" ] ; then
        awk 'BEGIN { OFS="\\t"} {print \$12, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' ${prefix}.genepred > ${prefix}.refflat
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '447'
    """
    touch ${prefix}.genepred
    touch ${prefix}.refflat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
