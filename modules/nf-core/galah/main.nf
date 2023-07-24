process GALAH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::galah=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/galah:0.3.1--hec16e2b_1':
        'biocontainers/galah:0.3.1--hec16e2b_1' }"

    input:
    tuple val(meta), path(bins)
    path(checkm_tab_table)
    path(genome_info)

    output:
    tuple val(meta), path("*.tsv")                    , emit: tsv
    tuple val(meta), path("${prefix}-dereplicated/*") , emit: dereplicated_bins, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def qc_input = checkm_tab_table ? "--checkm-tab-table ${checkm_tab_table}" : (genome_info ? "--genome-info ${genome_info}" : "")
    if (checkm_tab_table && genome_info) { error "genome_info table and checkm_tab_table both provided: please provide one or the other." }
    """
    mkdir ${prefix}-dereplicated

    galah cluster \\
        --threads ${task.cpus} \\
        --genome-fasta-files ${bins} \\
        ${qc_input} \\
        --output-cluster-definition ${prefix}-dereplicated_bins.tsv \\
        ${args}

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
