process PARABRICKS_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(interval_file)
    tuple val(ref_meta), path(fasta)

    output:
    tuple val(meta), path("*.vcf"),      emit: vcf,                 optional: true
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf,                optional: true
    path "compatible_versions.yml",      emit: compatible_versions, optional: true
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def output_file           = args.contains("--gvcf") ? "${prefix}.g.vcf.gz" : "${prefix}.vcf"
    def interval_file_command = interval_file ? interval_file.collect { "--interval-file ${it}" }.join(' ') : ""

    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun \\
        haplotypecaller \\
        --ref ${fasta} \\
        --in-bam ${input} \\
        --out-variants ${output_file} \\
        ${interval_file_command} \\
        ${num_gpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = args.contains("--gvcf") ? "echo '' | gzip > ${prefix}.g.vcf.gz" : "touch ${prefix}.vcf"
    """
    ${output_cmd}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun haplotypecaller --version 2>&1)

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
