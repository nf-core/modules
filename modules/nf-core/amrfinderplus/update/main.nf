process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.7--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' }"

    output:
    path "amrfinderdb.tar.gz", emit: db
    tuple val("${task.process}"), val('amrfinder'), eval('amrfinder --version'), emit: versions_amrfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    amrfinder_update -d amrfinderdb
    tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./
    """

    stub:
    """
    touch amrfinderdb.tar
    gzip amrfinderdb.tar
    """
}
