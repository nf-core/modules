process PARABRICKS_STARFUSION {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(genome_lib_dir)

    output:
    tuple val(meta), path(output_dir)   , emit: out_dir
    path "versions.yml"                 , emit: versions
    path "compatible_versions.yml"      , emit: compatible_versions, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        --ref \$INDEX \\
        --chimeric-junction ${chimeric_junction} \\
        --genome-lib-dir ${genome_lib_dir} \\
        --output-dir $prefix \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
        // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        --ref \$INDEX \\
        --chimeric-junction ${chimeric_junction} \\
        --genome-lib-dir ${genome_lib_dir} \\
        --output-dir $prefix \\
        $num_gpus \\
        $args

    # Capture once and build single-line compatible_with (spaces only, no tabs)
    pbrun_version_output=\$(pbrun fq2bam --version 2>&1)

    # Because of a space between BWA and mem in the version output this is handled different to the other modules
    compat_line=\$(echo "\$pbrun_version_output" | awk -F':' '
        /Compatible With:/ {on=1; next}
        /^---/ {on=0}
        on && /:/ {
            key=\$1; val=\$2
            gsub(/[ \\t]+/, " ", key); gsub(/^[ \\t]+|[ \\t]+\$/, "", key)
            gsub(/[ \\t]+/, " ", val); gsub(/^[ \\t]+|[ \\t]+\$/, "", val)
            a[++i]=key ": " val
        }
        END { for (j=1;j<=i;j++) printf "%s%s", (j>1?", ":""), a[j] }
    ')

    cat <<EOF > compatible_versions.yml
    "${task.process}":
    pbrun_version: \$(echo "\$pbrun_version_output" | awk '/^pbrun:/ {print \$2; exit}')
    compatible_with: "\$compat_line"
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
