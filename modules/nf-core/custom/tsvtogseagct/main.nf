process CUSTOM_TSVTOGSEAGCT {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.gct"), emit: gct

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    n_columns=\$(head -n 1 $tsv | tr "\t" "\n" | wc -l)
    n_lines=\$(wc -l < $tsv)
    gct_file=${prefix}.gct

    echo -e "#1.2\$(printf '\t%.0s' {1..\$n_columns})\n\$((n_lines-1))\t\$((n_columns-1))\$(printf '\t%.0s' {1..\$((n_columns-1))})" > \$gct_file
    echo -e "NAME\tDESCRIPTION\t\$(head -n 1 $tsv | cut -f1 --complement)" >> \$gct_file
    cat $tsv | tail -n +2 | awk 'BEGIN { OFS = "\t"} {\$1=\$1"\tNA"}1' >> \$gct_file
    """
}
