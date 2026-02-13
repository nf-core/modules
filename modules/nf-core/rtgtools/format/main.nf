process RTGTOOLS_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(input1), path(input2), path(sam_rg)

    output:
    tuple val(meta), path("*.sdf"), emit: sdf
    tuple val("${task.process}"), val('rgtools'), eval("rtg version | head -n 1 | sed 's/Product: RTG Tools //'"), topic: versions, emit: versions_rtgtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = meta.containsKey("single_end") ? meta.single_end : true

    def input = single ? "${input1}" : "--left ${input1} --right ${input2}"
    def rg = sam_rg ? "--sam-rg ${sam_rg}" : ""

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }

    """
    rtg RTG_MEM=${avail_mem} format \\
        ${args} \\
        ${rg} \\
        --output ${prefix}.sdf \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }
    """
    touch ${prefix}.sdf
    """
}
