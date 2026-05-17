process MIFASER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mifaser:1.64--pyh106432d_0' :
        'quay.io/biocontainers/mifaser:1.64--pyh106432d_0' }"

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*multi_ec.tsv"), emit: multi_ec
    tuple val(meta), path("*analysis.tsv"), emit: analysis
    tuple val(meta), path("*ec_count.tsv"), emit: ec_counts
    tuple val("${task.process}"), val('mi-faser'), eval("mifaser --version 2>&1 | sed 's/* v//'"), emit: versions_mifaser, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def input_flag = meta.single_end  ? "-f" : "-l"
    """
    mifaser \\
        ${args} \\
        ${input_flag} ${reads} \\
        --threads 1 \\
        --cpu ${task.cpus} \\
        --databasefolder \$PWD/${db} \\
        --outputfolder mifaser-${prefix}/

    for suf in multi_ec.tsv analysis.tsv ec_count.tsv; do
        mv mifaser-${prefix}/\${suf} ${prefix}_\${suf}
    done
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    for suf in multi_ec.tsv analysis.tsv ec_count.tsv; do
        touch ${prefix}_\${suf}
    done
    """
}
