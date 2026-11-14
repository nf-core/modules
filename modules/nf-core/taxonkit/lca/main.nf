// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAXONKIT_LCA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c33b5ce6cec42ed5986d11c706d5726b2eb7dc0902b9c1c333a319864a3b45ea/data'
        : 'community.wave.seqera.io/library/taxonkit:0.20.0--73c6bc587e9e62ac'}"

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
