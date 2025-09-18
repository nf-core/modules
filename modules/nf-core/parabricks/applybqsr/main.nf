process PARABRICKS_APPLYBQSR {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(bam_index)
    tuple val(meta3), path(bqsr_table)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "compatible_versions.yml", emit: compatible_versions, optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args             = task.ext.args    ?: ''
    def prefix           = task.ext.prefix  ?: "${meta.id}"
    def interval_command = intervals        ? intervals.collect { "--interval-file ${it}" }.join(' ') : ""
    def num_gpus         = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun \\
        applybqsr \\
        --ref ${fasta} \\
        --in-bam ${bam} \\
        --in-recal-file ${bqsr_table} \\
        ${interval_command} \\
        --out-bam ${prefix}.bam \\
        --num-threads ${task.cpus} \\
        ${num_gpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun applybqsr --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
