process TCOFFEE_EXTRACTFROMPDB {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4aed8f78ff751f9682f7f17c0259de2dafc57cee4ac5bd1aa0b2af8ba8ad495/data'
        : 'community.wave.seqera.io/library/t-coffee:13.46.0.919e8c6b--6242390d545250f1'}"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${prefix}.pdb"), emit: formatted_pdb
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TEMP='./'
    t_coffee -other_pg extract_from_pdb \
        -infile ${pdb} \
        $args \
        > "${prefix}.pdb"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch "${prefix}.pdb"
    """
}
