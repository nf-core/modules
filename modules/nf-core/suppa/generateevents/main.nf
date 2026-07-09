process SUPPA_GENERATEEVENTS {
    tag "${gtf}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    tuple val(meta), path(gtf)
    val format
    val pool_genes
    val event_type
    val boundary
    val threshold
    val exon_length

    output:
    tuple val(meta), path("*.{ioe,ioi,gtf}"), emit: events
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def pool_genes_arg = pool_genes ? "--pool-genes" : ''
    def event_type_arg = event_type ? "--event-type ${event_type}" : ''
    def boundary_arg = boundary ? "--boundary ${boundary}" : ''
    def threshold_arg = threshold ? "--threshold ${threshold}" : ''
    def exon_length_arg = exon_length ? "--exon-length ${exon_length}" : ''
    def ioe_args = format == 'ioe' ? [event_type_arg, boundary_arg, threshold_arg, exon_length_arg].join(' ') : ''

    """
    suppa.py \\
        generateEvents \\
        --input-file ${gtf} \\
        --format ${format} \\
        --output-file ${prefix} \\
        ${pool_genes_arg} \\
        ${ioe_args} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    echo ${args}
    touch ${prefix}.${format}
    """
}
