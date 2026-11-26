// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TREEMIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treemix:1.13--hffe338e_4' :
        'quay.io/biocontainers/treemix:1.13--hffe338e_4' }"

    input:
    tuple val(meta), path(input, arity: '1'), path(tree)

    output:
    tuple val(meta), path("*.llik")       , emit: llik
    tuple val(meta), path("*.treeout.gz") , emit: treeout
    tuple val(meta), path("*.vertices.gz"), emit: vertices
    tuple val(meta), path("*.edges.gz")   , emit: edges
    tuple val(meta), path("*.cov.gz")     , emit: cov
    tuple val(meta), path("*.modelcov.gz"), emit: modelcov
    tuple val(meta), path("*.covse.gz")   , optional: true, emit: covse
    tuple val(meta), path("*.log")        , emit: log
    tuple val("${task.process}"), val('treemix'), eval("treemix -h 2>&1 | sed -n 's/^TreeMix v[.] *//p'"), emit: versions_treemix, topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def tree_cmd = tree ? '-tf "' + tree + '"' : ''
    """
    treemix \\
        -i "${input}" \\
        -o "${prefix}" \\
        ${tree_cmd} \\
        ${args} \\
        > "${prefix}.log" 2>&1 \\
        || { cat "${prefix}.log" >&2; exit 1; }
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Exiting ln(likelihood) with 0 migration events: 0" > ${prefix}.llik
    echo "" | gzip > ${prefix}.treeout.gz
    echo "" | gzip > ${prefix}.vertices.gz
    echo "" | gzip > ${prefix}.edges.gz
    echo "" | gzip > ${prefix}.cov.gz
    echo "" | gzip > ${prefix}.modelcov.gz
    echo "" | gzip > ${prefix}.covse.gz
    printf 'TreeMix v. 1.13\\nSEED: 1\\n' > ${prefix}.log
    """
}
