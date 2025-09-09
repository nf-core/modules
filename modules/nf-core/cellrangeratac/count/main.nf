process CELLRANGERATAC_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger-atac:2.1.0"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}/outs/*"), emit: outs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGERATAC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def sample_arg = meta.samples.unique().join(",")
    def reference_name = reference.name
    """
    cellranger-atac \\
        count \\
        --id='${meta.id}' \\
        --fastqs=. \\
        --reference=$reference_name \\
        --sample=$sample_arg \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangeratac: \$(echo \$( cellranger-atac --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/{filtered_peak_bc_matrix.h5,filtered_tf_bc_matrix.h5,raw_peak_bc_matrix.h5,summary.csv,summary.json,peak_motif_mapping.bed,peak_annotation.tsv,cut_sites.bigwig,singlecell.csv}

    mkdir -p "${meta.id}/outs/filtered_peak_bc_matrix/"
    touch ${meta.id}/outs/filtered_peak_bc_matrix/{barcodes.tsv,peaks.bed,matrix.mtx}

    mkdir -p "${meta.id}/outs/filtered_tf_bc_matrix/"
    touch ${meta.id}/outs/filtered_tf_bc_matrix/motifs.tsv
    echo | gzip > "${meta.id}/outs/filtered_tf_bc_matrix/barcodes.tsv.gz"
    echo | gzip > "${meta.id}/outs/filtered_tf_bc_matrix/matrix.mtx.gz"

    mkdir -p "${meta.id}/outs/raw_peak_bc_matrix/"
    touch ${meta.id}/outs/raw_peak_bc_matrix/{barcodes.tsv,peaks.bed,matrix.mtx}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangeratac: \$(echo \$( cellranger-atac --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
