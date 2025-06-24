process CELLBENDER_REMOVEBACKGROUND {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'us.gcr.io/broad-dsde-methods/cellbender:0.3.2' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebcf140f995f79fcad5c17783622e000550ff6f171771f9fc4233484ee6f63cf/data':
        'community.wave.seqera.io/library/cellbender_webcolors:156d413fdfc16cdb' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5")               , emit: h5
    tuple val(meta), path("${prefix}_filtered.h5")      , emit: filtered_h5
    tuple val(meta), path("${prefix}_posterior.h5")     , emit: posterior_h5
    tuple val(meta), path("${prefix}_cell_barcodes.csv"), emit: barcodes
    tuple val(meta), path("${prefix}_metrics.csv")      , emit: metrics
    tuple val(meta), path("${prefix}_report.html")      , emit: report, optional: true
    tuple val(meta), path("${prefix}.pdf")              , emit: pdf
    tuple val(meta), path("${prefix}.log")              , emit: log
    tuple val(meta), path("ckpt.tar.gz")                , emit: checkpoint
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ""
    use_gpu = task.ext.use_gpu ? "--cuda" : ""
    """
    TMPDIR=. cellbender remove-background \
        ${args} \
        --cpu-threads ${task.cpus} \
        --estimator-multiple-cpu \
        ${use_gpu} \
        --input ${h5ad} \
        --output ${prefix}.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(cellbender --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5"
    touch "${prefix}_filtered.h5"
    touch "${prefix}_posterior.h5"
    touch "${prefix}_cell_barcodes.csv"
    touch "${prefix}_metrics.csv"
    touch "${prefix}_report.html"
    touch "${prefix}.pdf"
    touch "${prefix}.log"
    echo "" | gzip > ckpt.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(cellbender --version)
    END_VERSIONS
    """
}
