process CELLRANGER_AGGR {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    tuple val(meta), path(aggr_csv)

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_AGGR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    template "cellranger_aggr.py"

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_AGGR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    echo "$prefix" > ${prefix}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
