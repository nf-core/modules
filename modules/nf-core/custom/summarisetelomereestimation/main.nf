process CUSTOM_SUMMARISETELOMEREESTIMATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1fa34006114914735768188b781d7f2c8ae7132acd024ce4c45704170715b54f/data'
        : 'community.wave.seqera.io/library/pandas_python:8e99df08b7f7c3e1' }"

    input:
    tuple val(meta), path(length_tsv), path(content_tsv), val(length_tool)

    output:
    tuple val(meta), path("*_telomere_summary.tsv"), emit: summary
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'summarisetelomereestimation.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_telomere_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
