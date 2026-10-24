process GATK4_BEDTOINTERVALLIST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce8f4142326abbb74e28b97633a200b09439c54480803343b1cf13dc19e51ae5/data'
        : 'community.wave.seqera.io/library/gatk4-lite:4.7.0.0--79918b8b2632f4b1'}"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path('*.interval_list'), emit: interval_list
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK BedToIntervalList] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BedToIntervalList \\
        --INPUT ${bed} \\
        --OUTPUT ${prefix}.interval_list \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --TMP_DIR . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.interval_list
    """
}
