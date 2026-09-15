process SMUDGEPLOT_ALL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/320648514649f5379149eed196b162e0d409b785f670ddecdf948febd6917377/data':
        'community.wave.seqera.io/library/fastk_smudgeplot:1352fed7dbb39646' }"

    input:
    tuple val(meta), path(smu)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*_centralities.txt"), emit: centrality_txt
    tuple val(meta), path("*.smudge_report.tsv"), emit: smudge_report
    tuple val(meta), path("*.sma"), emit: sma
    tuple val("${task.process}"), val('smudgeplot'), eval('smudgeplot -v |& sed "s/.*v//"'), emit: versions_smudgeplot, topic: versions
    // FASTK does not report version to cli
    tuple val("${task.process}"), val('fastk'), val('1.2'), emit: versions_fastk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    // Export HOME to avoid issues with MATPLOTLIB needing a
    // writable config directory
    """
    export HOME=\$PWD/nxf_home

    smudgeplot all \\
        ${args} \\
        -o ${prefix} \\
        ${smu}
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.smu
    """
}
