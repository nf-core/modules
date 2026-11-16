// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GAPPA_EXAMINEASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a30e4444b69c99470cae4a631a9e293ed76d475a6a44812a1f38ad2b9ece7af/data'
        : 'community.wave.seqera.io/library/gappa:0.9.0--2f0efa59bcb34481'}"

    input:
    tuple val(meta), path(jplace), path(taxonomy)

    output:
    tuple val(meta), path("*profile.tsv")         , emit: profile
    tuple val(meta), path("*labelled_tree.newick"), emit: labelled_tree
    tuple val(meta), path("*per_query.tsv")       , emit: per_query, optional: true
    tuple val(meta), path("*krona.profile")       , emit: krona    , optional: true
    tuple val(meta), path("*sativa.tsv")          , emit: sativa   , optional: true
    tuple val("${task.process}"), val('gappa'), eval("gappa --version 2>&1 | sed 's/v//'"), emit: versions_gappa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // gappa reads a gzipped jplace natively, but a gzipped taxon file is read as text and
    // fails with "A line in the taxon file didn't have two tab separated columns".
    def taxonfile = taxonomy.name.endsWith('.gz') ? taxonomy.baseName : "${taxonomy}"
    def gunzip    = taxonomy.name.endsWith('.gz') ? "gzip -cd ${taxonomy} > ${taxonfile}" : ""
    """
    $gunzip

    gappa \\
        examine assign \\
        ${args} \\
        --threads ${task.cpus} \\
        --jplace-path ${jplace} \\
        --taxon-file ${taxonfile} \\
        --file-prefix ${prefix}.
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.tsv
    touch ${prefix}.labelled_tree.newick
    """
}
