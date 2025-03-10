process CUSTOM_TABULARTOGSEACLS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(samples)

    output:
    tuple val(meta), path("*.cls"), emit: cls
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(mawk -W version | head -n 1 | awk '{print  \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cls
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(mawk -W version | head -n 1 | awk '{print  \$2}')
    END_VERSIONS
    """

}
