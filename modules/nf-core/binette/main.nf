process BINETTE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/binette:1.1.2--pyh7e72e81_0'
        : 'biocontainers/binette:1.1.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(input_binning, stageAs: "input_binning/*"), path(contigs), path(proteins)
    val input_type
    path checkm2_db

    output:
    tuple val(meta), path("${meta.id}_final_bins/*")                  , emit: refined_bins
    tuple val(meta), path("${meta.id}_final_bins_quality_reports.tsv"), emit: refined_bins_report
    tuple val(meta), path("${meta.id}_input_bins_quality_reports/*")  , emit: input_bins_report
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def protein_arg = proteins ? "--proteins ${proteins}" : ""
    def input_arg   = ""

    if (input_type == 'fasta') {
        input_arg = "--bin_dirs"
    } else if (input_type == 'tsv') {
        input_arg = "--contig2bin_tables"
    } else {
        error "Invalid input_type: ${input_type}. Must be 'fasta' or 'tsv'"
    }

    """
    binette \\
        ${input_arg} input_binning/* \\
        ${protein_arg} \\
        --contigs ${contigs} \\
        --threads ${task.cpus} \\
        --checkm2_db ${checkm2_db} \\
        -o . \\
        ${args}

    mv final_bins "${prefix}_final_bins"
    mv final_bins_quality_reports.tsv ${prefix}_final_bins_quality_reports.tsv
    mv input_bins_quality_reports ${prefix}_input_bins_quality_reports

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}_final_bins
    touch ${prefix}_final_bins/bin_b.fasta

    mkdir ${prefix}_input_bins_quality_reports
    touch ${prefix}_input_bins_quality_reports/bin_a.tsv
    touch ${prefix}_final_bins_quality_reports.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """
}
