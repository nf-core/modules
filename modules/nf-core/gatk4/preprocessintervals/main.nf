process GATK4_PREPROCESSINTERVALS {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1330cb37b3acde07db0a04befa76973fea72d1d10a79542132e3d77950225af/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:b945148cc5e2a616'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(dict)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(exclude_intervals)

    output:
    tuple val(meta), path("${prefix}.interval_list"), emit: interval_list
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def include_command = intervals ? "--intervals ${intervals}" : ""
    def exclude_command = exclude_intervals ? "--exclude-intervals ${exclude_intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PreprocessIntervals] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PreprocessIntervals \\
        ${include_command} \\
        ${exclude_command} \\
        --reference ${fasta} \\
        --output ${prefix}.interval_list \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.interval_list
    """
}
