process OPT_STAT {
    tag "$meta.id"
    label 'process_high'

    container "khersameesh24/opt:v0.0.1"

    input:
    tuple val(meta), path(probe_targets)
    tuple val(meta2), path(fwd_oriented_probes)
    path(gene_synonyms)

    output:
    tuple val(meta), path("${meta.id}/collapsed_summary.tsv"), emit: summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "OPT_STAT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def synonyms = gene_synonyms ? "-s ${gene_synonyms}": ""

    """
    opt \\
    -o ${prefix} \\
    stat \\
    -i ${probe_targets} \\
    -q ${fwd_oriented_probes} \\
    ${synonyms} \\
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
    touch "${prefix}/collapsed_summary.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opt: \$(opt --version)
    END_VERSIONS
    """
}
