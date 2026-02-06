process PARABRICKS_INDEXGVCF {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path("*.{idx,tbi}"), emit: gvcf_index
    path "compatible_versions.yml",       emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun \\
        indexgvcf \\
        --input ${gvcf} \\
        ${num_gpus} \\
        ${args}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = gvcf.any { item -> item.name.endsWith(".gz") } ? "touch ${prefix}.g.vcf.gz.tbi" : "touch ${prefix}.g.vcf.idx"
    """
    ${output_cmd}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun indexgvcf --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
