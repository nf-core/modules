process CUSTOM_TABULARTOGSEACHIP {

    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(tabular)
    tuple val(id)  , val(symbol)

    output:
    tuple val(meta), path("*.chip"), emit: chip
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    function find_column_number {
        file=\$1
        column=\$2

        head -n 1 \$file | tr '\\t' '\\n' | grep -n "^\${column}\$" | awk -F':' '{print \$1}'
    }

    id_col=\$(find_column_number $tabular $id)
    symbol_col=\$(find_column_number $tabular $symbol)
    outfile=\$(echo $tabular | sed 's/\\(.*\\)\\..*/\\1/').chip

    echo -e "Probe Set ID\\tGene Symbol\\tGene Title" > \${outfile}.tmp
    tail -n +2 $tabular | awk -F'\\t' -v id=\$id_col -v symbol=\$symbol_col '{print \$id"\\t"\$symbol"\\tNA"}' >> \${outfile}.tmp
    mv \${outfile}.tmp \${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """

    stub:
    """
    outfile=\$(echo $tabular | sed 's/\\(.*\\)\\..*/\\1/').chip
    touch \$outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
