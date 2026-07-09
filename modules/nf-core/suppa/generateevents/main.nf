process SUPPA_GENERATEEVENTS {
    tag "${gtf}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    path gtf
    val format
    val pool_genes
    val event_type
    val boundary
    val threshold
    val exon_length

    output:
    path "events*.*", emit: events
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def pg = pool_genes ? "--pool-genes" : ''
    def et = event_type ? "--event-type ${event_type}" : ''
    def bd = boundary ? "--boundary ${boundary}" : ''
    def th = threshold ? "--threshold ${threshold}" : ''
    def el = exon_length ? "--exon-length ${exon_length}" : ''
    def ioe_args = [et, bd, th, el].join(' ')

    if (format == 'ioe') {
    """
    suppa.py \\
        generateEvents \\
        --input-file ${gtf} \\
        --format ${format} \\
        --output-file events \\
        ${pg} \\
        ${ioe_args} \\
        ${args}
    """
    } else if (format == 'ioi') {
    """
    suppa.py \\
        generateEvents \\
        --input-file ${gtf} \\
        --format ${format} \\
        --output-file events \\
        ${pg} \\
        ${args}
    """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    echo ${args}

    if [[ "${format}" == "ioe" ]]; then
        if [[ "${pool_genes}" == "true" ]]; then
            touch "events.ioe"
        else
            IFS=' ' read -ra events <<< "${event_type}"
            for event in "\${events[@]}"; do
                touch "events_\${event}.ioe"
            done
        fi
    elif [[ "${format}" == "ioi" ]]; then
        touch "events.ioi"
    fi
    """
}
