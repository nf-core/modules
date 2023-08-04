process GALAH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::galah=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/galah%3A0.3.1--h031d066_3':
        'biocontainers/galah:0.3.1--h031d066_3' }"

    input:
    tuple val(meta), path(bins), path(qc_table), val(qc_format)

    output:
    tuple val(meta), path("*.tsv")                    , emit: tsv
    tuple val(meta), path("${prefix}-dereplicated/*") , emit: dereplicated_bins
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def qc_args = (qc_format == "checkm") ? "--checkm-tab-table ${qc_table}" : "--genome-info ${qc_table}"
    def qc_input = qc_table ? qc_args : ""
    def valid_qc_format = qc_format in ["checkm", "genome_info"]
    if( qc_table && !valid_qc_format ) {
        error "Invalid qc_format supplied! qc_format should be either 'checkm' or 'genome_info'."
    }
    """
    mkdir ${prefix}-dereplicated

    galah cluster \\
        --threads ${task.cpus} \\
        --genome-fasta-files ${bins} \\
        ${qc_input} \\
        --output-cluster-definition ${prefix}-dereplicated_bins.tsv \\
        --output-representative-fasta-directory ${prefix}-dereplicated

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version | sed 's/galah //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}-dereplicated/
    touch ${prefix}-dereplicated/test.fa
    touch ${prefix}-dereplicated_bins.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version | sed 's/galah //')
    END_VERSIONS
    """
}
