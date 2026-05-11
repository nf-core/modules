process MSISENSORPRO_BASELINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro%3A1.3.0--hfef96ef_0':
        'biocontainers/msisensor-pro:1.3.0--hfef96ef_0' }"

    input:
    tuple val(meta), path(list)
    tuple val(meta2), path(pro_all_files, stageAs: 'pro_all/*')

    output:
    tuple val(meta), path("${prefix}.baseline"), emit: baseline
    tuple val("${task.process}"), val('msisensor-pro'), eval("msisensor-pro --version 2>&1 | sed -nE 's/Version:\\s*//p'") , emit: versions_msisensorpro, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build the configure file consumed by \`msisensor-pro baseline -i\`:
    #   one row per normal sample, columns: <sample_name>\\t<absolute_path_to_*_all>
    # Sample names are derived from each file's basename with the "_all" suffix stripped.
    : > configure.txt
    for f in pro_all/*_all; do
        [ -e "\$f" ] || continue
        name=\$(basename "\$f" _all)
        printf '%s\\t%s\\n' "\$name" "\$(readlink -f "\$f")" >> configure.txt
    done

    msisensor-pro \\
        baseline \\
        -d $list \\
        -i configure.txt \\
        -o ${prefix}.baseline \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.baseline
    """
}
