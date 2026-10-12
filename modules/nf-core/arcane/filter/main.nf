process ARCANE_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcane:1.0.0--pyh106432d_0':
        'quay.io/biocontainers/arcane:1.0.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*/*_genes.txt"),                                 emit: genes
    tuple val(meta), path("*/*_mapping_gene_ids.pickle"),                   emit: pickle
    tuple val(meta), path("*/*_new_id_to_original.csv"),                    emit: csv
    tuple val(meta), path("*/*_ref.fa.gz"),                                 emit: reference
    tuple val(meta), path("*/*_extracted_genes.fa.gz"),                     emit: extracted_genes, optional: true
    tuple val("${task.process}"), val('arcane'), eval("arcane --version"),  topic: versions, emit: versions_arcane

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    arcane filter \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        --prefix ${prefix} \\
        --name ${prefix} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo "" | gzip > ${prefix}/arcane_${prefix}_ref.fa.gz
    echo "" | gzip > ${prefix}/arcane_${prefix}_extracted_genes.fa.gz
    touch ${prefix}/arcane_${prefix}_genes.txt
    touch ${prefix}/arcane_${prefix}_new_id_to_original.csv
    touch ${prefix}/arcane_${prefix}_mapping_gene_ids.pickle
    """
}
