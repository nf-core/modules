process PARABRICKS_DEEPVARIANT {
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
    tuple val(meta), path(input), path(input_index), path(interval_file)
    path  fasta

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def copy_index_command = input_index ? "cp -L $input_index `readlink -f $input`.bai" : ""
    """
    $copy_index_command

    pbrun \\
        deepvariant \\
        --ref $fasta \\
        --in-bam $input \\
        --out-variants ${prefix}.vcf \\
        $interval_file_command \\
        --num-gpus $task.accelerator.request \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}

// TODO
// gvcf output options 
// additional deepvariant arguments
