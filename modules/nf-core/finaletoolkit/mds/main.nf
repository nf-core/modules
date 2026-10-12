process FINALETOOLKIT_MDS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34f01d128ed135aedc33ddb62fced3911bef6d1a909694291b7184bf83719402/data'
        : 'community.wave.seqera.io/library/finaletoolkit:0.11.1--8fe5ba6ec9e2ec95'}"

    input:
    tuple val(meta), path(motifs)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val("${task.process}"), val('finaletoolkit'), eval("finaletoolkit --version | sed 's/FinaleToolkit //g'"), topic: versions, emit: versions_finaletoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_mds"

    // NOTE: finaletoolkit outputs to stdout without any other information
    // To make it more useful, this module redirects the outputs to a TSV

    """

    export MPLCONFIGDIR=.

    echo -e "sample\tmds" > ${prefix}.tsv

    finaletoolkit \\
        mds \\
        ${motifs} \\
        ${args} | \\
        sed 's/^/${meta.id}\t/1' >> ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_mds"
    """
    echo ${args}

    touch "${prefix}.tsv"
    """
}
