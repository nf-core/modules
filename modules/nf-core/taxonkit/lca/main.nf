process TAXONKIT_LCA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.20.0--h9ee0642_0':
        'quay.io/biocontainers/taxonkit:0.20.0--h9ee0642_0' }"

    input:
    tuple val(meta), val(taxids), path(taxidsfile)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('taxonkit'), eval("taxonkit version | sed 's/.* v//'"), emit: versions_taxonkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!taxids && !taxidsfile) {
            error("Error: Neither taxids nor taxidsfile provided! ")
    }    
    if (taxids && taxidsfile) {
            error("Error: Both taxids and taxidsfile provided, and only one should be passed! ")
    }       
    """
    taxonkit \\
        lca \\
        ${args} \\
        --data-dir ${taxdb} \\
        --threads ${task.cpus} \\
        --out-file ${prefix}.tsv \\
        ${taxids? "<<< '$taxids'": taxidsfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
