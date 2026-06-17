process GT_LTRHARVEST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'quay.io/biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*.tabout")       , emit: tabout      , optional: true
    tuple val(meta), path("*.gff3")         , emit: gff3        , optional: true
    tuple val(meta), path("$out_name")      , emit: fasta       , optional: true    // When args has -out
    tuple val(meta), path("$outinner_name") , emit: inner_fasta , optional: true    // When args has -outinner
    tuple val("${task.process}"), val('genometools'), eval("gt --version | sed '1!d;s/.* //'"), emit: versions_gt, topic: versions


    script:
    def args        = task.ext.args                 ?: ''
    def prefix      = task.ext.prefix               ?: "${meta.id}"
    def extension   = args.contains("-tabout no")   ? "gff3" : "tabout"
    out_name        = (args.split('-').find { arg_it -> arg_it =~ /out .*\.(fa|fsa|fasta)/ } ?: 'out.fasta').replace('out ', '').trim()
    outinner_name   = (args.split('-').find { arg_it -> arg_it =~ /outinner .*\.(fa|fsa|fasta)/ } ?: 'outinner.fasta').replace('outinner ', '').trim()
    """
    gt \\
        ltrharvest \\
        -index ${index}/suffixerator \\
        ${args} \\
        > ${prefix}.${extension}
    """

    stub:
    def args        = task.ext.args                 ?: ''
    def prefix      = task.ext.prefix               ?: "${meta.id}"
    def extension   = args.contains("-tabout no")   ? "gff3"                        : "tabout"

    out_name        = (args.split('-').find { arg_it -> arg_it =~ /out .*\.(fa|fsa|fasta)/ } ?: 'out.fasta').replace('out ', '').trim()
    outinner_name   = (args.split('-').find { arg_it -> arg_it =~ /outinner .*\.(fa|fsa|fasta)/ } ?: 'outinner.fasta').replace('outinner ', '').trim()

    def touch_out   = args.contains("-out")         ? "touch $out_name"             : ''
    def touch_inner = args.contains("-outinner")    ? "touch $outinner_name"        : ''
    """
    touch "${prefix}.${extension}"
    ${touch_out}
    ${touch_inner}
    """
}
