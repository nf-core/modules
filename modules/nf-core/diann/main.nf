process DIANN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/diann/v1.8.1_cv1/diann_v1.8.1_cv1.img' :
        'docker.io/biocontainers/diann:v1.8.1_cv1' }"

    input:
    tuple val(meta), path(ms_files), val(ms_file_names), path(fasta), path(library), path(quant, stageAs: 'quant/*')

    output:
    // Library outputs
    tuple val(meta), path("${prefix}.predicted.speclib"), emit: predict_speclib, optional: true
    tuple val(meta), path("${prefix}.speclib"), emit: final_speclib, optional: true
    tuple val(meta), path("${prefix}.tsv.skyline.speclib"), emit: skyline_speclib, optional: true
    
    // Quantification outputs
    tuple val(meta), path("*.quant"), emit: diann_quant, optional: true
    
    // Report outputs (from final quantification)
    tuple val(meta), path("${prefix}.tsv"), emit: main_report, optional: true
    tuple val(meta), path("${prefix}.parquet"), emit: report_parquet, optional: true
    tuple val(meta), path("${prefix}.manifest.txt"), emit: report_manifest, optional: true
    tuple val(meta), path("${prefix}.protein_description.tsv"), emit: protein_description, optional: true
    tuple val(meta), path("${prefix}.stats.tsv"), emit: report_stats, optional: true
    tuple val(meta), path("${prefix}.pr_matrix.tsv"), emit: pr_matrix, optional: true
    tuple val(meta), path("${prefix}.pg_matrix.tsv"), emit: pg_matrix, optional: true
    tuple val(meta), path("${prefix}.gg_matrix.tsv"), emit: gg_matrix, optional: true
    tuple val(meta), path("${prefix}.unique_genes_matrix.tsv"), emit: unique_gene_matrix, optional: true
    
    // Common outputs
    tuple val(meta), path("*.log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}" ?: "diann"
    
    // Handle MS files input: two modes depending on whether we need actual files or just names
    // - ms_files: Actual file paths - used when DIA-NN needs to read raw MS data
    // - ms_file_names: Just basenames - used with --use-quant when DIA-NN only needs file names
    //   to match against preprocessed .quant files in quant/ directory, avoiding unnecessary file staging
    def ms_input = ''
    if (ms_files && ms_files != []) {
        ms_input = ms_files instanceof List ? ms_files.collect{ "--f ${it}" }.join(' ') : "--f ${ms_files}"
    } else if (ms_file_names && ms_file_names != []) {
        ms_input = ms_file_names instanceof List ? ms_file_names.collect{ "--f ${it}" }.join(' ') : "--f ${ms_file_names}"
    }
    
    def fasta_input = fasta && fasta != [] ? "--fasta ${fasta}" : ''
    def lib_input = library && library != [] ? "--lib ${library}" : ''

    // When quant files are provided, set temp directory and enable --use-quant
    // These flags must be used together: --temp points to quant files, --use-quant tells DIA-NN to use them
    def quant_args = quant && quant != [] ? "--temp ./quant/ --use-quant" : "--temp ./"

    """
    diann \\
        ${ms_input} \\
        ${fasta_input} \\
        ${lib_input} \\
        --threads ${task.cpus} \\
        --out-lib ${prefix} \\
        --out ${prefix}.tsv \\
        ${quant_args} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}" ?: "diann"

    """
    # Library outputs
    touch ${prefix}.predicted.speclib
    touch ${prefix}.speclib
    touch ${prefix}.tsv
    touch ${prefix}.tsv.skyline.speclib
    
    # Quant outputs
    touch ${prefix}.quant
    
    # Report outputs
    touch ${prefix}.tsv
    touch ${prefix}.parquet
    touch ${prefix}.manifest.txt
    touch ${prefix}.protein_description.tsv
    touch ${prefix}.stats.tsv
    touch ${prefix}.pr_matrix.tsv
    touch ${prefix}.pg_matrix.tsv
    touch ${prefix}.gg_matrix.tsv
    touch ${prefix}.unique_genes_matrix.tsv
    
    # Common outputs
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: 1.8.1
    END_VERSIONS
    """
}
