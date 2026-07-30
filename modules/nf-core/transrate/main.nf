process TRANSRATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transrate:1.0.3--h87e0c26_7':
        'quay.io/biocontainers/transrate:1.0.3--h87e0c26_7' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_mqc.csv")                                              , emit: multiqc
    tuple val("${task.process}"), val('transrate'), eval("transrate --version"), topic: versions, emit: versions_transrate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    transrate \\
        $args \\
        --threads $task.cpus \\
        --assembly $assembly \\
        --output ${prefix}_transrate

    mv ${prefix}_transrate/assemblies.csv ${prefix}_mqc.csv

    # transrate writes the resolved absolute path of the input assembly into the first
    # column of the report; replace it with the sample prefix so the report content is
    # reproducible across environments and task work directories.
    sed -i "2s@^[^,]*@${prefix}@" ${prefix}_mqc.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.csv
    """
}
