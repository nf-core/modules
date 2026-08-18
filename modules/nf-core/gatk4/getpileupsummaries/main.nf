process GATK4_GETPILEUPSUMMARIES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1330cb37b3acde07db0a04befa76973fea72d1d10a79542132e3d77950225af/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:b945148cc5e2a616'}"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path variants
    path variants_tbi

    output:
    tuple val(meta), path('*.pileups.table'), emit: table
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals ${intervals}" : "--intervals ${variants}"
    def reference_command = fasta ? "--reference ${fasta}" : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GetPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GetPileupSummaries \\
        --input ${input} \\
        --variant ${variants} \\
        --output ${prefix}.pileups.table \\
        ${reference_command} \\
        ${interval_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileups.table
    """
}
