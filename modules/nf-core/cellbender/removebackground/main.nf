process CELLBENDER_REMOVEBACKGROUND {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cellbender:0.3.0--c4addb97ab2d83fe':
        'community.wave.seqera.io/library/cellbender:0.3.0--41318a055fc3aacb' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5")               , emit: h5
    tuple val(meta), path("${prefix}_filtered.h5")      , emit: filtered_h5
    tuple val(meta), path("${prefix}_posterior.h5")     , emit: posterior_h5
    tuple val(meta), path("${prefix}_cell_barcodes.csv"), emit: barcodes
    tuple val(meta), path("${prefix}_metrics.csv")      , emit: metrics
    tuple val(meta), path("${prefix}_report.html")      , emit: report
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
    touch "ckpt.tar.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(cellbender --version)
    END_VERSIONS
    """
}
