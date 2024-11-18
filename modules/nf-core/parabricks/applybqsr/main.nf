process PARABRICKS_APPLYBQSR {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.2-1"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bam_index)
    tuple val(meta3), path(bqsr_table)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? intervals.collect{"--interval-file $it"}.join(' ') : ""
    def copy_index_command = bam_index ? "cp -L $bam_index `readlink -f $bam`.bai" : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    # parabricks complains when index is not a regular file in the same directory as the bam
    # copy the index to this path.
    $copy_index_command

    pbrun \\
        applybqsr \\
        --ref $fasta \\
        --in-bam $bam \\
        --in-recal-file $bqsr_table \\
        $interval_command \\
        --out-bam ${prefix}.bam \\
        --num-threads $task.cpus \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? intervals.collect{"--interval-file $it"}.join(' ') : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
