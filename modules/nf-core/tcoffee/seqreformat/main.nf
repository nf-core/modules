process TCOFFEE_SEQREFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.46.0.919e8c6b--hfc96bf3_0':
        'quay.io/biocontainers/t-coffee:13.46.0.919e8c6b--hfc96bf3_0' }"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: formatted_file
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TEMP='./'
    t_coffee -other_pg seq_reformat \
        -in ${infile} \
        $args \
        > "${prefix}.txt"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch "${prefix}.txt"
    """
}
