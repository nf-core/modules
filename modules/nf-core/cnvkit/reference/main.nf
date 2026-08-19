process CNVKIT_REFERENCE {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/3935f8d7507f85fd20e073f0625bb6bdff8aa1b6044f4bbb720c9a8063ea0390/data'
:         'community.wave.seqera.io/library/cnvkit:0.9.14--288e98d6210b7304' }"

    input:
    path fasta
    path targets
    path antitargets

    output:
    path "*.cnn",        emit: cnn
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    cnvkit.py \\
        reference \\
        --fasta ${fasta} \\
        --targets ${targets} \\
        --antitargets ${antitargets} \\
        --output ${prefix}.reference.cnn \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    touch ${prefix}.reference.cnn
    """
}
