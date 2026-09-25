process VCONTACT3_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcontact3:3.1.6--pyhdfd78af_1':
        'quay.io/biocontainers/vcontact3:3.1.6--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(genome)
    path database

    output:
    tuple val(meta), path("${prefix}/"), emit: results
    tuple val("${task.process}"), val('vcontact3'), eval('vcontact3 version'), emit: versions_vcontact3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    vcontact3 run \\
        -n ${genome} \\
        -o ${prefix}/ \\
        --threads ${task.cpus} \\
        -d ${database} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/clusters.csv
    touch ${prefix}/merged_df.csv
    touch ${prefix}/vContact_contig_id.txt
    """
}
