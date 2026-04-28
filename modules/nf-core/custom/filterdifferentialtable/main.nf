process CUSTOM_FILTERDIFFERENTIALTABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b06950ac325030db5976f3d9c536e358eb686503af35c7e6222f86d016b3466f/data' :
        'community.wave.seqera.io/library/pandas_python:67bda66f0cb8a241' }"

    input:
    tuple val(meta), path(input_file)
    tuple val(logfc_column), val(fc_threshold), val(fc_cardinality)
    tuple val(stat_column), val(stat_threshold), val(stat_cardinality)

    output:
    tuple val(meta), path("*_filtered.tsv")     , emit: filtered
    tuple val(meta), path("*_filtered_up.tsv")  , emit: filtered_up
    tuple val(meta), path("*_filtered_down.tsv"), emit: filtered_down
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python

    import pandas as pd
    import sys
    from math import log2

    if not "${input_file}".lower().endswith(('.csv', '.tsv', '.txt')):
        sys.exit("Please provide a .csv, .tsv or .txt file!")

    # Determine the separator based on file extension
    sep = ',' if "${input_file}".lower().endswith('.csv') else '\t'

    # Read the input file
    table = pd.read_csv("${input_file}", sep=sep)

    # Calculate log2 fold change threshold
    logfc_threshold = log2(float("${fc_threshold}"))

    # define evaluation
    def evaluate_condition(x, threshold, cardinality):
        if cardinality == ">=":
            return x >= threshold
        elif cardinality == "<=":
            return x <= threshold
        elif cardinality == ">":
            return x > threshold
        elif cardinality == "<":
            return x < threshold
        else:
            raise ValueError(f"Invalid cardinality: {cardinality}")

    # Apply filters
    mask = (
        table["${logfc_column}"].notna() &
        table["${stat_column}"].notna() &
        table["${logfc_column}"].abs().apply(lambda x: evaluate_condition(x, logfc_threshold, "${fc_cardinality}")) &
        table["${stat_column}"].apply(lambda x: evaluate_condition(x, float("${stat_threshold}"), "${stat_cardinality}"))
    )
    filtered_table = table[mask]
    filtered_table_up = filtered_table[filtered_table["${logfc_column}"] > 0]
    filtered_table_down = filtered_table[filtered_table["${logfc_column}"] < 0]

    # Write the filtered table
    filtered_table.to_csv("${prefix}_filtered.tsv", sep='\t', index=False)
    filtered_table_up.to_csv("${prefix}_filtered_up.tsv", sep='\t', index=False)
    filtered_table_down.to_csv("${prefix}_filtered_down.tsv", sep='\t', index=False)

    # Write versions
    with open('versions.yml', 'w') as version_file:
        version_file.write('"${task.process}":\\n')
        version_file.write(f"    pandas: {pd.__version__}\\n")
        version_file.write(f"    python: {sys.version.split()[0]}\\n")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.tsv
    touch ${prefix}_filtered_up.tsv
    touch ${prefix}_filtered_down.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(pip show pandas | sed -n 's/Version: //p')
        python: \$(python --version | cut -f 2 -d " ")
    END_VERSIONS
    """
}
