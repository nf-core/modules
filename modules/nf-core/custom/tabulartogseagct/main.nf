process CUSTOM_TABULARTOGSEAGCT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

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
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
