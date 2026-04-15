process PARSNP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/parsnp:2.1.5--0605933fc69e7b20':
        'community.wave.seqera.io/library/parsnp:2.1.5--2c7f64ad14a79523' }"

    input:
    tuple val(meta), path(genomes, stageAs: "genomes")
    path reference

    output:
    tuple val(meta), path("*.xmfa")        , emit: xmfa
    tuple val(meta), path("*.ggr")         , emit: ggr
    tuple val(meta), path("*.snps.mblocks"), emit: snps_mblocks, optional: true
    tuple val(meta), path("*.tree")        , emit: tree
    tuple val(meta), path("partition")     , emit: partition, optional: true
    tuple val("${task.process}"), val('parsnp'), eval("parsnp --version 2>/dev/null | tail -n 1 | sed 's/parsnp //'"), topic: versions, emit: versions_parsnp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parsnp \\
        -r "${reference}" \\
        -d "${genomes}" \\
        -o parsnp_out \\
        -p ${task.cpus} \\
        ${args}

    mv parsnp_out/parsnp.xmfa ${prefix}.xmfa
    mv parsnp_out/parsnp.ggr ${prefix}.ggr
    if [ -f parsnp_out/parsnp.snps.mblocks ]; then
        mv parsnp_out/parsnp.snps.mblocks "${prefix}.snps.mblocks"
    fi
    mv parsnp_out/parsnp.tree ${prefix}.tree
    if [ -d parsnp_out/partition ]; then
        mv parsnp_out/partition .
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xmfa
    touch ${prefix}.ggr
    touch ${prefix}.snps.mblocks
    touch ${prefix}.tree
    """
}
