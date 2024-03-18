process CELLRANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:7.1.0"

    input:
    tuple val(meta), path(reads, stageAs: "fastq_???/*")
    path  reference

    output:
    tuple val(meta), path("**/outs/**")                          , emit: outs
    tuple val(meta), path("**/outs/filtered_feature_bc_matrix**"), emit: filtered
    tuple val(meta), path("**/outs/raw_feature_bc_matrix**")     , emit: raw
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template "cellranger_count.py"

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    mkdir -p "${prefix}/outs/filtered_feature_bc_matrix"
    mkdir -p "${prefix}/outs/raw_feature_bc_matrix"
    echo "$prefix" > ${prefix}/outs/fake_file.txt
    echo "$prefix" > ${prefix}/outs/filtered_feature_bc_matrix/fake_file.txt
    echo "$prefix" > ${prefix}/outs/raw_feature_bc_matrix/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
