process GALAH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/galah:0.4.2--h7b50bb2_1':
        'biocontainers/galah:0.4.2--h7b50bb2_1' }"

    input:
    tuple val(meta), path(bins), path(qc_table), val(qc_format)

    output:
    tuple val(meta), path("*.tsv")      , emit: tsv
    tuple val(meta), path("${prefix}/*"), emit: dereplicated_bins
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def qc_input = ""
    if(qc_format == "checkm2") {
        qc_input = "--checkm2-quality-report ${qc_table}"
    } else if(qc_format == "checkm") {
        qc_input = "--checkm-tab-table ${qc_table}"
    } else if(qc_format == "genome-info") {
        qc_input = "--genome-info ${qc_table}"
    }
    """
    mkdir -p ${prefix}-dereplicated

    galah cluster \\
        --threads ${task.cpus} \\
        --genome-fasta-files ${bins} \\
        ${qc_input} \\
        --output-cluster-definition ${prefix}.tsv \\
        --output-representative-fasta-directory ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version | sed 's/galah //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}/
    touch ${prefix}/test.fa
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version | sed 's/galah //')
    END_VERSIONS
    """
}
