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

    output:
    tuple val(meta), path("${prefix}"), emit: db
    tuple val("${task.process}"), val("centrifuger"), eval("centrifuger -v 2>&1 | head -n 1 | cut -d ' ' -f 2"), emit: versions_centrifuger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // Always use -l option: input FASTA file(s) are converted into a reference list
    def refs = reference instanceof List ? reference : [reference]
    // Create an input str so we can be used to create a file
    def reference_list_str = refs.join("\n")

    // check if conversion table is given.
    if (!conversion_table) {
        error "centrifuger-build requires a --conversion-table"
    }

   """
    mkdir -p ${prefix}

    #Create reference file from input files --> use it with -l option
    cat <<EOF > reference_list.txt
    ${reference_list_str}
    EOF

    centrifuger-build \\
        -l reference_list.txt \\
        --taxonomy-tree $taxonomy_nodes \\
        --name-table $taxonomy_names \\
        --conversion-table ${conversion_table} \\
        -t ${task.cpus} \\
        -o ${prefix}/${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir -p ${prefix}/
    touch ${prefix}/${prefix}.1.cfr
    touch ${prefix}/${prefix}.2.cfr
    touch ${prefix}/${prefix}.3.cfr
    touch ${prefix}/${prefix}.4.cfr
    """
}
