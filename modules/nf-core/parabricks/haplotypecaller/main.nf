process PARABRICKS_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"

    /*
    NOTE: Parabricks requires the files to be non-symlinked
    Do not change the stageInMode to soft linked! This is default on Nextflow.
    If you change this setting be careful.
    */
    stageInMode "copy"

    input:
    tuple val(meta), path(bam), path(bam_index), path(intervals)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

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

    """
    # parabricks complains when index is not a regular file in the same directory as the bam
    # copy the index to this path.
    $copy_index_command

    pbrun \\
        haplotypecaller \\
        --ref $fasta \\
        --in-bam $bam \\
        $interval_command \\
        --out-variants ${prefix}.vcf.gz
        --num-threads $task.cpus \\
        --num-gpus $task.accellerator.request \\
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
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
