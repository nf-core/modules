// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process UNIVERSC {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/nf-core/universc:1.2.7"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val(technology)

    output:
    tuple val(meta), path("${prefix}/outs/*"), emit: outs
    tuple val("${task.process}"), val('cellranger'), eval("cellranger --version | sed 's/cellranger //'"), emit: versions_cellranger, topic: versions
    tuple val("${task.process}"), val('universc'), eval("sed -n 's/^universcversion=\"\\(.*\\)\"/\\1/p' \$(which launch_universc.sh)"), emit: versions_universc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args        = task.ext.args   ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--file $reads" : "-R1 ${reads[0]} -R2 ${reads[1]}"

    def reference_name = reference.name
    """
    # UniverSC edits its own and Cell Ranger's install directories, so run from task-local copies
    cr=\$(dirname \$(which cellranger))
    mkdir -p .cellranger/lib .cellranger/mro
    cp -r \$cr/lib/python .cellranger/lib/
    cp -r \$cr/mro/rna .cellranger/mro/
    printf '#!/bin/sh\\nexec %s/cellranger "\$@"\\n' \$cr > .cellranger/cellranger
    chmod +x .cellranger/cellranger
    cp -r \$(dirname \$(which launch_universc.sh)) .universc
    export PATH="\$PWD/.cellranger:\$PWD/.universc:\$PATH"

    export PYTHON_EGG_CACHE=\$(pwd)/.cache
    launch_universc.sh \\
        --id ${prefix} \\
        ${input_reads} \\
        --technology ${technology} \\
        --reference ${reference_name} \\
        --jobmode "local" \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --per-cell-data \\
        ${args}

    # save log files
    echo !! > ${prefix}/outs/_invocation
    """


    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "UNIVERSC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/outs/
    cd ${prefix}/outs/

    touch _invocation

    touch basic_stats.txt
    touch metrics_summary.csv
    touch molecule_info.h5
    touch possorted_genome_bam.bam
    touch possorted_genome_bam.bam.bai
    touch web_summary.html

    mkdir -p filtered_feature_bc_matrix
    touch filtered_feature_bc_matrix.h5
    echo "" | gzip > filtered_feature_bc_matrix/barcodes.tsv.gz
    echo "" | gzip > filtered_feature_bc_matrix/features.tsv.gz
    echo "" | gzip > filtered_feature_bc_matrix/matrix.mtx.gz

    mkdir -p raw_feature_bc_matrix
    touch raw_feature_bc_matrix.h5
    echo "" | gzip > raw_feature_bc_matrix/barcodes.tsv.gz
    echo "" | gzip > raw_feature_bc_matrix/features.tsv.gz
    echo "" | gzip > raw_feature_bc_matrix/matrix.mtx.gz

    cd ../..
    """
}
