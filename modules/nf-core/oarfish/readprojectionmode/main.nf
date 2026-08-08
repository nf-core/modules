process OARFISH_READPROJECTIONMODE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/890d90ebb61b3756527fabe481d801a8e203453fde05a67ccc01640a9a445cd9/data':
        'community.wave.seqera.io/library/oarfish:0.10.3--837f9667a3c87dbb' }"

    input:
    tuple val(meta), path(reads), val(seq_tech)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(annotation)

    output:
    tuple val(meta), path("${prefix}.quant")         , emit: quant
    tuple val(meta), path("${prefix}.meta_info.json"), emit: meta_info
    tuple val(meta), path("${prefix}.ambig_info.tsv"), emit: ambig_info
    tuple val(meta), path("${prefix}.infreps.pq")    , optional: true, emit: infreps
    tuple val(meta), path("${prefix}.prob.lz4")      , optional: true, emit: prob
    tuple val("${task.process}"), val('oarfish'), eval("oarfish --version | sed 's/oarfish //'"), topic: versions, emit: versions_oarfish

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    oarfish \\
        --reads ${reads} \\
        --genome ${fasta} \\
        --annotation ${annotation} \\
        --seq-tech ${seq_tech} \\
        -o ${prefix} \\
        --threads ${task.cpus} \\
        --filter-group no-filters \\
        --verbose \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.quant
    touch ${prefix}.ambig_info.tsv
    touch ${prefix}.infreps.pq
    touch ${prefix}.prob.lz4
    touch ${prefix}.meta_info.json
    """
}
