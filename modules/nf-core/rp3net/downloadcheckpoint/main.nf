process RP3NET_DOWNLOADCHECKPOINT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.4' :
        'quay.io/biocontainers/wget:1.21.4' }"

    output:
    path("*.ckpt"), emit: checkpoint
    tuple val("${task.process}"), val('wget'), eval("wget --version 2>&1 | head -1 | sed 's/GNU Wget //' | sed 's/ .*//'"), topic: versions, emit: versions_wget

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    wget \\
        -q \\
        ${args} \\
        https://ftp.ebi.ac.uk/pub/software/RP3Net/v0.1/checkpoints/rp3net_v0.1_d.ckpt
    """

    stub:
    """
    touch rp3net_v0.1_d.ckpt
    """
}
