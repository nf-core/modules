process ARCANE_EXPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcane:1.0.0--pyh106432d_0':
        'quay.io/biocontainers/arcane:1.0.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(r1), path(r2)
    tuple val(meta2), path(hash), path(info)
    val chemistry

    output:
    tuple val(meta), path("*_counts.mex.gz"),                               emit: counts
    tuple val(meta), path("*_genes.tsv.gz"),                                emit: genes
    tuple val(meta), path("*_barcodes.tsv.gz"),                             emit: barcodes
    tuple val("${task.process}"), val('arcane'), eval("arcane --version"),  emit: versions_arcane, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_base = hash.baseName
    """
    arcane express \
        --index ${index_base} \
        --R1 ${r1} \
        --R2 ${r2} \
        --chemistry ${chemistry} \
        --out ${prefix} \
        --threads ${task.cpus} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_counts.mex.gz
    echo "" | gzip > ${prefix}_genes.tsv.gz
    echo "" | gzip > ${prefix}_barcodes.tsv.gz
    """
}
