process GATK4_SHIFTFASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*_shift.fasta"), emit: shift_fa
    tuple val(meta), path("*_shift.fasta.fai"), emit: shift_fai
    tuple val(meta), path("*_shift.back_chain"), emit: shift_back_chain
    tuple val(meta), path("*_shift.dict"), emit: dict, optional: true
    tuple val(meta), path("*.intervals"), emit: intervals, optional: true
    tuple val(meta), path("*.shifted.intervals"), emit: shift_intervals, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_dict = dict ? "--sequence-dictionary ${dict}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ShiftFasta] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ShiftFasta \\
        --reference ${fasta} \\
        --output ${prefix}_shift.fasta \\
        --shift-back-output ${prefix}_shift.back_chain \\
        ${args} \\
        ${seq_dict} \\
        --tmp-dir .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.intervals
    touch ${prefix}_shift.back_chain
    touch ${prefix}_shift.dict
    touch ${prefix}.shifted.intervals
    touch ${prefix}_shift.fasta
    touch ${prefix}_shift.fasta.fai
    """
}
