process PARABRICKS_RNAFQ2BAM {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(index)
    tuple val(meta3), path(genome_lib_dir)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq ${reads}" : "--in-fq ${reads}"
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp ${fasta} \$INDEX

    pbrun \\
        rna_fq2bam  \\
        --ref \$INDEX \\
        ${in_fq_command} \\
        --output-dir ${prefix} \\
        --genome-lib-dir ${genome_lib_dir} \\
        --out-bam ${prefix}.bam \\
        ${num_gpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    ${qc_metrics_output}
    ${duplicate_metrics_output}

    # Capture the full version output once and store it in a variable
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
