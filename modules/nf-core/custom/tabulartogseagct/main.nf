// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_TABULARTOGSEAGCT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d809369ef0dfb6f580ebb812d33e442db2866d73757acc6a490ad9aa1a55ef10/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--e94c8172b8bba3cf'}"

    input:
    tuple val(meta), path(tabular)

    output:
    tuple val(meta), path("*.gct"), emit: gct
    tuple val("${task.process}"), val('mawk'), eval("mawk --version | sed '1!d;s/mawk //;s/ .*//'"), emit: versions_mawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: []
    def prefix = task.ext.prefix ?: "${meta.id}"
    def separator = args.separator ? "${args.separator}" : (tabular.getName().endsWith(".csv") ? ',' : '\t')
    separator = separator == '\t' ? '\\t' : separator

    """
    n_columns=\$(head -n 1 ${tabular} | tr "${separator}" "\\n" | wc -l)
    n_lines=\$(wc -l < ${tabular})
    gct_file=${prefix}.gct

    echo -e "#1.2\$(printf '\\t%.0s' {1..\$n_columns})\\n\$((n_lines-1))\\t\$((n_columns-1))\$(printf '\\t%.0s' {1..\$((n_columns-1))})" > \$gct_file
    echo -e "NAME\\tDESCRIPTION\\t\$(head -n 1 ${tabular} | cut -f1 -d\$'${separator}' --complement | awk -F'${separator}' 'BEGIN { OFS = "\\t"}; {\$1=\$1}1' )" >> \$gct_file
    cat ${tabular} | tail -n +2 | awk -F'${separator}' 'BEGIN { OFS = "\\t"} {\$1=\$1"\\tNA"}1' >> \$gct_file
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gct
    """
}
