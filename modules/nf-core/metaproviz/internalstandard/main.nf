process METAPROVIZ_INTERNALSTANDARD {
    tag "$meta.id"
    label 'process_single'

    // Shared MetaProViz image (built from modules/nf-core/metaproviz/Dockerfile).
    container "ghcr.io/saezlab/metaproviz:0.0.1"

    input:
    // Two acceptable, mutually exclusive shapes: se_rds alone (a
    // SummarizedExperiment .rds), or all three of data_matrix/
    // feature_matrix/sample_matrix (plain TSV/CSV). Pass [] for whichever
    // shape you're not using — see tests/main.nf.test for both call shapes.
    tuple val(meta), path(se_rds), path(data_matrix), path(feature_matrix), path(sample_matrix)

    output:
    // Read-only QC: reports internal-standard CV, does NOT modify the SE.
    // every output now carries meta (except versions.yml,
    // which never does in nf-core convention) so downstream channels keep
    // the sample id attached to every file, not just cv_table.
    tuple val(meta), path("*.cv.tsv")            , emit: cv_table
    tuple val(meta), path("*.high_var.txt")      , emit: high_var
    tuple val(meta), path("*.condition_cv.tsv")  , emit: condition_cv
    tuple val(meta), path("*.plots.rds")         , emit: plots
    tuple val(meta), path("*.report.html")       , emit: report
    tuple val(meta), path("*.log")               , emit: log
    path "versions.yml"                                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def se_arg  = se_rds ? "--se ${se_rds}" : ""
    def csv_arg = (data_matrix && feature_matrix && sample_matrix) ?
        "--data_matrix ${data_matrix} --feature_matrix ${feature_matrix} --sample_matrix ${sample_matrix}" : ""
    """
    internal_standard.R \\
        $se_arg \\
        $csv_arg \\
        $args \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
        metaproviz: \$(Rscript -e 'cat(as.character(packageVersion("MetaProViz")))')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.cv.tsv ${prefix}.high_var.txt \\
          ${prefix}.condition_cv.tsv ${prefix}.plots.rds \\
          ${prefix}.report.html ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.1
        metaproviz: $VERSION
    END_VERSIONS
    """
}
