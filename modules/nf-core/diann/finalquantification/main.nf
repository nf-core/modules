process DIANN_FINALQUANTIFICATION {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/diann/v1.8.1_cv1/diann_v1.8.1_cv1.img' :
        'docker.io/biocontainers/diann:v1.8.1_cv1' }"

    input:
    tuple val(meta), val(ms_files), path(fasta), path(empirical_library), path("quant/")

    output:
    tuple val(meta), path("diann_report.tsv"), emit: main_report, optional: true
    tuple val(meta), path("diann_report.parquet"), emit: report_parquet, optional: true
    tuple val(meta), path("diann_report.manifest.txt"), emit: report_manifest, optional: true
    tuple val(meta), path("diann_report.protein_description.tsv"), emit: protein_description, optional: true
    tuple val(meta), path("diann_report.stats.tsv"), emit: report_stats
    tuple val(meta), path("diann_report.pr_matrix.tsv"), emit: pr_matrix
    tuple val(meta), path("diann_report.pg_matrix.tsv"), emit: pg_matrix
    tuple val(meta), path("diann_report.gg_matrix.tsv"), emit: gg_matrix
    tuple val(meta), path("diann_report.unique_genes_matrix.tsv"), emit: unique_gene_matrix
    tuple val(meta), path("diannsummary.log"), emit: log
    tuple val(meta), path("empirical_library.tsv"), emit: final_speclib, optional: true
    tuple val(meta), path("empirical_library.tsv.skyline.speclib"), emit: skyline_speclib, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_FINALQUANTIFICATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''

    """
    # Notes: if .quant files are passed, mzml/.d files are not accessed, so the name needs to be passed but files
    # do not need to pe present.

    diann   --lib ${empirical_library} \\
            --fasta ${fasta} \\
            --f ${(ms_files as List).join(' --f ')} \\
            --threads ${task.cpus} \\
            --temp ./quant/ \\
            ${args}

    cp diann_report.log.txt diannsummary.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_FINALQUANTIFICATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''

    """
    touch diann_report.tsv
    touch diann_report.parquet
    touch diann_report.manifest.txt
    touch diann_report.protein_description.tsv
    touch diann_report.stats.tsv
    touch diann_report.pr_matrix.tsv
    touch diann_report.pg_matrix.tsv
    touch diann_report.gg_matrix.tsv
    touch diann_report.unique_genes_matrix.tsv
    touch diannsummary.log
    touch empirical_library.tsv
    touch empirical_library.tsv.skyline.speclib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """
}
