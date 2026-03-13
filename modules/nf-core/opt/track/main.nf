process OPT_TRACK {
    tag "$meta.id"
    label 'process_high'

    container "khersameesh24/opt:v0.0.1"

    input:
    tuple val(meta), path(fwd_oriented_fa)
    tuple val(meta2), path(ref_annot_gff), path(ref_annot_fa)

    output:
    tuple val(meta), path("${meta.id}/probe2targets.tsv"), emit: probes2target
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "OPT_TRACK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    opt \\
    -o ${prefix} \\
    -p ${task.cpus} \\
    track \\
    -q ${fwd_oriented_fa} \\
    -a ${ref_annot_gff} \\
    -t ${ref_annot_fa} \\
    ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opt: \$(opt --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/probe2targets.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opt: \$(opt --version)
    END_VERSIONS
    """
}
