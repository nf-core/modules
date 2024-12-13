process CUSTOM_TABULARTOGSEAGCT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/coreutils:8.30--b947da103164f84b"

    input:
    tuple val(meta), path(tabular)

    output:
    tuple val(meta), path("*.gct"), emit: gct
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: []
    def prefix = task.ext.prefix ?: "${meta.id}"
    def separator = args.separator ? "${args.separator}" : ( tabular.getName().endsWith(".csv") ? ',': '\t' )
    separator = separator == '\t' ? '\\t': separator

    """
    n_columns=\$(head -n 1 $tabular | tr "$separator" "\\n" | wc -l)
    n_lines=\$(wc -l < $tabular)
    gct_file=${prefix}.gct

    echo -e "#1.2\$(printf '\\t%.0s' {1..\$n_columns})\\n\$((n_lines-1))\\t\$((n_columns-1))\$(printf '\\t%.0s' {1..\$((n_columns-1))})" > \$gct_file
    echo -e "NAME\\tDESCRIPTION\\t\$(head -n 1 $tabular | cut -f1 -d\$'$separator' --complement | awk -F'$separator' 'BEGIN { OFS = "\\t"}; {\$1=\$1}1' )" >> \$gct_file
    cat $tabular | tail -n +2 | awk -F'$separator' 'BEGIN { OFS = "\\t"} {\$1=\$1"\\tNA"}1' >> \$gct_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(mawk -W version | head -n 1 | awk '{print  \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gct
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(mawk -W version | head -n 1 | awk '{print  \$2}')
    END_VERSIONS
    """
}
