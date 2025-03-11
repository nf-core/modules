process GATK4_SHIFTFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*_shift.fasta")       , emit: shift_fa
    tuple val(meta), path("*_shift.fasta.fai")   , emit: shift_fai
    tuple val(meta), path("*_shift.back_chain")  , emit: shift_back_chain
    tuple val(meta), path("*_shift.dict")        , emit: dict              , optional: true
    tuple val(meta), path("*.intervals")         , emit: intervals         , optional: true
    tuple val(meta), path("*.shifted.intervals") , emit: shift_intervals   , optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_dict = dict ? "--sequence-dictionary ${dict}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ShiftFasta] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ShiftFasta \\
        --reference $fasta \\
        --output ${prefix}_shift.fasta \\
        --shift-back-output ${prefix}_shift.back_chain \\
        $args \\
        $seq_dict \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch test.intervals
    touch test_shift.back_chain
    touch test_shift.dict
    touch test.shifted.intervals
    touch test_shift.fasta
    touch test_shift.fasta.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
