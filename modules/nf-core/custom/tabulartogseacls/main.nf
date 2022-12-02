process CUSTOM_TABULARTOGSEACLS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(samples)

    output:
    tuple val(meta), path("*.cls"), emit: cls
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def separator = samples.getName().endsWith(".csv") ? ',' : '\\t'
    """
    cls_file=${prefix}.cls

    column_number=\$(cat $samples | head -n 1 | tr '$separator' "\\n" | grep -En "^$meta.variable" | awk -F':' '{print \$1}')
    classes=\$(tail -n +2 $samples | awk -F'$separator' '{print \$'\$column_number'}')
    unique_classes=\$(echo -e "\$classes" | awk '!x[\$0]++')

    echo -e "\$(echo -e \"\$classes\" | wc -l) \$(echo -e \"\$unique_classes\" | wc -l) 1" > \$cls_file
    echo -e "#\$(echo -e \"\$unique_classes\" | tr '\\n' ' ')" | sed "s/ \$//" >> \$cls_file
    echo -e "\$classes" | tr '\\n' ' ' | sed "s/ \$//" >> \$cls_file
    echo -e "\\n" >> \$cls_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
