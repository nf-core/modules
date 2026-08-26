process FINALETOOLKIT_FRAGLENGTHBINS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34f01d128ed135aedc33ddb62fced3911bef6d1a909694291b7184bf83719402/data':
        'community.wave.seqera.io/library/finaletoolkit:0.11.1--8fe5ba6ec9e2ec95' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.frag_len.tsv"), emit: tsv
    tuple val(meta), path("*_hist.png"), emit: png
    tuple val("${task.process}"), val('finaletoolkit'), eval("finaletoolkit --version | sed 's/FinaleToolkit //g'"), topic: versions, emit: versions_finaletoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    finaletoolkit \\
        frag-length-bins \\
        ${bam} \\
        ${args} \\
        --histogram-path ${prefix}_hist.png \\
        -o "${prefix}.frag_len.tsv"
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.frag_len.tsv
    touch ${prefix}_hist.png
    """
}
