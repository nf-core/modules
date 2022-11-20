process PARABRICKS_APPLYBQSR {
    tag "$meta.id"
    label 'process_high'
    memory "30 GB"
    cpus 8
    accelerator 1

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used with Parabricks at the moment. Please use a Docker container."
    }
    if (workflow.containerEngine == 'singularity') {
        exit 1, "Singularity containers cannot be used with Parabricks at the moment. Please use a Docker container."
    }

    // container "nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1"
    container "645946264134.dkr.ecr.us-west-2.amazonaws.com/clara-parabricks:4.0.0-1"

    input:
    tuple val(meta), path(bam), path(fasta), path(recal_file)
    path interval_file

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    """
    pbrun \\
        applybqsr \\
        --ref $fasta \\
        --in-bam $bam \\
        --in-recal-file $recal_file \\
        $interval_file_command \\
        --out-bam ${prefix}_bqsr.bam \\
        --num-threads $task.cpus \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

