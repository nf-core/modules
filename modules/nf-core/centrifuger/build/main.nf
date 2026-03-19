process CENTRIFUGER_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuger:1.1.0--hf426362_0':
        'biocontainers/centrifuger:1.1.0--hf426362_0' }"
    input:
    tuple val(meta), path(reference)
    path taxonomy_nodes
    path taxonomy_names
    path conversion_table
    path reference_list



    output:
    tuple val(meta), path("${prefix}/"), emit: db
    tuple val("${task.process}"), val("centrifuger"), eval("centrifuger -v 2>&1 | head -n 1 | cut -d ' ' -f 2"), emit: versions_centrifuger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // if reference list inserted and not conversiont able --> use -l
    def reference_input = (reference_list && !conversion_table)
        ? "-l ${reference_list}"
        : reference.collect { r -> "-r ${r}" }.join(' ')
    def conversion_arg = conversion_table ? "--conversion-table ${conversion_table}" : ''
   """
    mkdir -p ${prefix}

    centrifuger-build \\
        ${reference_input} \\
        --taxonomy-tree $taxonomy_nodes \\
        --name-table $taxonomy_names \\
        ${conversion_arg} \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/${prefix}.1.cfr
    touch ${prefix}/${prefix}.2.cfr
    touch ${prefix}/${prefix}.3.cfr
    touch ${prefix}/${prefix}.4.cfr
    """
}
