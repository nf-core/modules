// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process EPANG_FILTERQUERY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/0600fdb9287986d3d684b95d6e94f597c41b24a746608532d7c0836bd7b44c15/data'
:         'community.wave.seqera.io/library/gawk_gzip:fe0c202000bc8ed4' }"

    input:
    tuple val(meta), path(reference, stageAs: 'reference/*'), path(query, stageAs: 'query/*'), val(min_query_columns), val(min_reference_coverage)

    output:
    tuple val(meta), path("*.query.fasta.gz"), emit: query, optional: true
    tuple val(meta), path("*.filterquery.tsv"), emit: table
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    assert min_query_columns as Integer >= 1 : "EPANG_FILTERQUERY: min_query_columns must be at least 1, the lowest value EPA-ng accepts (got ${min_query_columns})"
    assert (min_reference_coverage as BigDecimal) >= 0 && (min_reference_coverage as BigDecimal) <= 1 : "EPANG_FILTERQUERY: min_reference_coverage must be between 0 and 1 (got ${min_reference_coverage})"
    template 'filterquery.sh'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.query.fasta.gz
    printf 'query\\treference_columns\\tquery_columns\\treference_coverage\\tkept\\n' > ${prefix}.filterquery.tsv
    """
}
