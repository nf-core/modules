process CELLRANGER_AGGR {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    tuple val(meta), val(sample_ids)
    path(molecule_h5_files, stageAs: "sample_???/*")

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    tuple val("${task.process}"), val('cellranger'), eval("cellranger --version 2>&1 | sed 's/.*cellranger-//'"), emit: versions_cellranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_AGGR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_ids_str = sample_ids.join(' ')
    """
    echo "sample_id,molecule_h5" > aggregation.csv
    sids=(${sample_ids_str})
    dirs=(\$(ls -d sample_* | sort))
    for i in "\${!sids[@]}"; do
        h5=\$(find "\${dirs[\$i]}" -name "*molecule_info.h5" | head -1)
        echo "\${sids[\$i]},\${h5}" >> aggregation.csv
    done

    cellranger aggr \\
        --id "${prefix}" \\
        --csv aggregation.csv \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_AGGR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    echo "$prefix" > ${prefix}/outs/fake_file.txt
    touch aggregation.csv
    """
}
