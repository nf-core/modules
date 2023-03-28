process CELLRANGER_ARC_COUNT {
    tag "$meta.id"
    label 'process_high'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_COUNT module does not support Conda.
        Please use docker or singularity containers."
    }
    container "heylf/cellranger-arc:2.0.2"

    input:
    tuple val(meta), path(reads)
    path  lib_csv
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs/*"), emit: outs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def lib_csv_name = lib_csv.name
    def reference_name = reference.name
    """
    cellranger-arc \\
        count \\
        --id='${meta.id}' \\
        --libraries=$lib_csv_name \\
        --reference=$reference_name \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
