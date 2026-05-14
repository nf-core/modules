process HUMANN3_REGROUP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.6.1--pyh7cba7a3_0' :
        'quay.io/biocontainers/humann:3.6.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(input)
    val groups
    path utility_db

    output:
    tuple val(meta), path("*_regroup.tsv.gz"), emit: regroup
    tuple val("${task.process}"), val('HUMAnN'), eval("humann --version 2>&1 | sed 's/humann v//'"), emit: versions_humann, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ $input == *.gz ]]; then
        gunzip -c $input > input.tsv
    else
        cp $input input.tsv
    fi

    printf '%s\n' \\
        'import os, sys' \\
        'import humann.config as config' \\
        'config.utility_mapping_database = os.environ["HUMANN_UTILITY_DB"]' \\
        'from humann.tools.regroup_table import main' \\
        'sys.exit(main())' \\
        > run_regroup.py

    export HUMANN_UTILITY_DB="${utility_db}"

    python run_regroup.py \\
        --input input.tsv \\
        --output ${prefix}_regroup.tsv \\
        --groups ${groups} \\
        ${args}

    gzip -n ${prefix}_regroup.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_regroup.tsv.gz
    """
}
