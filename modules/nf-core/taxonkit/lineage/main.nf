// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAXONKIT_LINEAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47a901198dfa6fd00c447557e953a64abb8794d49812ccb623bc973510069706/data'
        : 'community.wave.seqera.io/library/taxonkit:0.18.0--36a4cff30b2fab4e'}"

    input:
    tuple val(meta), val(taxid), path(taxidfile)
    path taxdb

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('taxonkit'), eval("taxonkit version | sed 's/.* v//'"), emit: versions_taxonkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert (!taxid && taxidfile) || (taxid && !taxidfile)
    """
    taxonkit \\
        lineage \\
        ${args} \\
        --data-dir ${taxdb} \\
        --threads ${task.cpus} \\
        --out-file ${prefix}.tsv \\
        ${taxid? "<<< '${taxid}'": taxidfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
