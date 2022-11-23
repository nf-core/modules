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

    container "nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1"
    // container "645946264134.dkr.ecr.us-west-2.amazonaws.com/clara-parabricks:4.0.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(interval_file)
    path  fasta

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def copy_index_command = input_index ? "cp -L $input_index `readlink -f $input`.bai" : ""
    """
    # parabricks complains when index is not a regular file in the same directory as the bam
    # copy the index to this path. 
    $copy_index_command

    pbrun \\
        applybqsr \\
        --ref $fasta \\
        --in-bam $input \\
        --in-recal-file $bqsr_table \\
        $interval_file_command \\
        --out-bam ${prefix}.bam \\
        --num-threads $task.cpus \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

// TODO: limit to 2 gpus (applybqsr won't use more than 2, according to the documentation)
