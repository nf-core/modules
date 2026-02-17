process BEDOPS_GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.42--h9948957_0':
        'biocontainers/bedops:2.4.42--h9948957_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    tuple val("${task.process}"), val("bedops"), eval('bedops --version | sed -n "s/.*version: *\\([^ ]*\\).*/\\1/p"'), emit: versions_bedops, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    cat \\
    $gtf \\
    | gtf2bed \\
    $args \\
    --attribute-key=exon_id \\
    > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.bed
    """

}
