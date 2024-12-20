process CUSTOM_FILTERDIFFERENTIALTABLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(input_file)
    val(logFC_column)
    val(FC_threshold)
    val(padj_column)
    val(padj_threshold)

    output:
    tuple val(meta), path("*_filtered.tsv"), emit: filtered
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: []
    def prefix = task.ext.prefix ?: "${meta.id}"
    def comparison1 = args.comparison1 ?: ">="
    def comparison2 = args.comparison2 ?: "<="
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
    logFC_threshold = log2(float("${FC_threshold}"))

    # define evaluation
    def evaluate_condition(x, threshold, comparison):
        if comparison == ">=":
            return x >= threshold
        elif comparison == "<=":
            return x <= threshold
        else:
            raise ValueError(f"Invalid comparison: {comparison}")

    # Apply filters
    mask = (
        table["${logFC_column}"].notna() &
        table["${padj_column}"].notna() &
        table["${logFC_column}"].abs().apply(lambda x: evaluate_condition(x, logFC_threshold, "${comparison1}")) &
        table["${padj_column}"].apply(lambda x: evaluate_condition(x, float("${padj_threshold}"), "${comparison2}"))
    )
    filtered_table = table[mask]

    # Write the filtered table
    filtered_table.to_csv("${prefix}_filtered.tsv", sep='\t', index=False)

    # Write versions
    with open('versions.yml', 'w') as version_file:
        version_file.write('"${task.process}":\\n')
        version_file.write(f"    pandas: {pd.__version__}\\n")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered.tsv
    echo '"${task.process}":\\n    pandas: 1.5.2' > versions.yml
    """
}
