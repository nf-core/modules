process PARABRICKS_MINIMAP2 {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(interval_file)
    tuple val(meta4), path(known_sites)
    val output_fmt

    output:
    tuple val(meta), path("*.bam"),                   emit: bam,                 optional: true
    tuple val(meta), path("*.bai"),                   emit: bai,                 optional: true
    tuple val(meta), path("*.cram"),                  emit: cram,                optional: true
    tuple val(meta), path("*.crai"),                  emit: crai,                optional: true
    tuple val(meta), path("*.table"),                 emit: bqsr_table,          optional: true
    tuple val(meta), path("*_qc_metrics"),            emit: qc_metrics,          optional: true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics,   optional: true
    path "compatible_versions.yml",                   emit: compatible_versions, optional: true
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Determine input format and set appropriate command flag
    def in_command = ""
    if (reads.name.matches(".*\\.(fastq|fq)(\\.gz)?\$")) {
        in_command = "--in-fq ${reads}"
    }
    else if (reads.name.endsWith('.bam')) {
        in_command = "--in-bam ${reads}"
    }
    else {
        error("Unsupported input file format: ${reads.name}. Supported formats: fastq, fastq.gz, fq, fq.gz, bam")
    }
    def extension = "${output_fmt}"

    def known_sites_command    = known_sites ? (known_sites instanceof List ? known_sites.collect { "--knownSites ${it}" }.join(' ') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command  = interval_file ? (interval_file instanceof List ? interval_file.collect { "--interval-file ${it}" }.join(' ') : "--interval-file ${interval_file}") : ""

    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun \\
        minimap2 \\
        --ref ${fasta} \\
        ${in_command} \\
        --out-bam ${prefix}.${extension} \\
        ${known_sites_command} \\
        ${known_sites_output_cmd} \\
        ${interval_file_command} \\
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
    def extension       = "${output_fmt}"
    def extension_index = "${output_fmt}" == "cram" ? "crai" : "bai"
    def known_sites_output       = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output        = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.${extension}
    touch ${prefix}.${extension}.${extension_index}
    ${known_sites_output}
    ${qc_metrics_output}
    ${duplicate_metrics_output}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun minimap2 --version 2>&1)

    # We handle this different to the other modules because minimap does not begin with an Uppercase letter

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | tr '\\t' ' ' | awk -F':' '/Compatible With:/,/^---/ { if (\$0 !~ /Compatible With:/ && \$0 !~ /^---\$/ && index(\$0,":")>0) { key=\$1; val=\$2; gsub(/^[ ]+|[ ]+\$/, "", key); gsub(/^[ ]+|[ ]+\$/, "", val); printf "  %s: %s\\n", key, val } }')
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
