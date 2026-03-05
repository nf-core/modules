process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7150110bea918bd523e41aab58f3be2eda949bbe7f715500681b01965e5f667/data':
        'community.wave.seqera.io/library/gfastats:1.3.11--692f8595aa6b06a8' }"

    input:
    tuple val(meta), path(assembly)
    val out_fmt
    val genome_size
    val target
    tuple val(meta2), path(agpfile)
    tuple val(meta3), path(include_bed)
    tuple val(meta4), path(exclude_bed)
    tuple val(meta5), path(instructions)

    output:
    tuple val(meta), path("*.assembly_summary"), emit: assembly_summary
    tuple val(meta), path("*.${out_fmt}.gz")   , emit: assembly        , optional: true
    tuple val("${task.process}"), val('gfastats'), eval("gfastats -v | sed '1!d;s/.*v//'"), emit: versions_gfastats, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agp    = agpfile ? "--agp-to-path $agpfile" : ""
    def ibed   = include_bed ? "--include-bed $include_bed" : ""
    def ebed   = exclude_bed ? "--exclude-bed $exclude_bed" : ""
    def sak    = instructions ? "--swiss-army-knife $instructions" : ""
    def output_sequences = out_fmt ? "--out-format ${prefix}.${out_fmt}.gz" : ""
    """
    gfastats \\
        $args \\
        --threads $task.cpus \\
        $agp \\
        $ibed \\
        $ebed \\
        $sak \\
        $output_sequences \\
        --input-sequence $assembly \\
        $genome_size \\
        $target \\
        > ${prefix}.assembly_summary
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.${out_fmt}.gz
    touch ${prefix}.assembly_summary
    """
}
