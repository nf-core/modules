process CUSTOM_TABULARTOGSEACLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d809369ef0dfb6f580ebb812d33e442db2866d73757acc6a490ad9aa1a55ef10/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--e94c8172b8bba3cf'}"

    input:
    tuple val(meta), path(samples)

    output:
    tuple val(meta), path("*.cls"), emit: cls
    tuple val("${task.process}"), val('mawk'), eval("mawk --version | sed '1!d;s/mawk //;s/ .*//'"), emit: versions_mawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: []
    def prefix = task.ext.prefix ?: "${meta.id}"
    def separator = args.separator ? "${args.separator}" : (samples.getName().endsWith(".tsv") ? '\t' : ',')
    separator = separator == '\t' ? '\\t' : separator
    def variable = args.variable
    if (!variable) {
        error("Supply a variable in the sample sheet from which to derive classes")
    }
    """
    cls_file=${prefix}.cls

    column_number=\$(cat ${samples} | head -n 1 | tr '${separator}' "\\n" | grep -En "^${variable}\$" | awk -F':' '{print \$1}')
    classes=\$(tail -n +2 ${samples} | awk -F'${separator}' '{print \$'\$column_number'}' | sed 's/^\$/empty/g')
    unique_classes=\$(echo -e "\$classes" | awk '!x[\$0]++')

    echo -e "\$(echo -e \"\$classes\" | wc -l) \$(echo -e \"\$unique_classes\" | wc -l) 1" > \$cls_file
    echo -e "#\$(echo -e \"\$unique_classes\" | tr '\\n' ' ')" | sed "s/ \$//" >> \$cls_file
    echo -e "\$classes" | tr '\\n' ' ' | sed "s/ \$//" >> \$cls_file
    echo -e "\\n" >> \$cls_file
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cls
    """

}
