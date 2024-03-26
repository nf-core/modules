process GT_LTRHARVEST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.5--py310h3db02ab_0':
        'biocontainers/genometools-genometools:1.6.5--py310h3db02ab_0' }"

    input:
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*.tabout")       , emit: tabout      , optional: true
    tuple val(meta), path("*.gff3")         , emit: gff3        , optional: true
    tuple val(meta), path("$out_name")      , emit: fasta       , optional: true    // When args has -out
    tuple val(meta), path("$outinner_name") , emit: inner_fasta , optional: true    // When args has -outinner
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args                 ?: ''
    def prefix      = task.ext.prefix               ?: "${meta.id}"
    def extension   = args.contains("-tabout no")   ? "gff3" : "tabout"
    out_name        = (args.split('-').find { it =~ /out .*\.(fa|fsa|fasta)/ } ?: 'out.fasta').replace('out ', '').trim()
    outinner_name   = (args.split('-').find { it =~ /outinner .*\.(fa|fsa|fasta)/ } ?: 'outinner.fasta').replace('outinner ', '').trim()
    """
    gt \\
        ltrharvest \\
        -index "$index/suffixerator" \\
        $args \\
        > "${prefix}.${extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args                 ?: ''
    def prefix      = task.ext.prefix               ?: "${meta.id}"
    def extension   = args.contains("-tabout no")   ? "gff3"                        : "tabout"

    out_name        = (args.split('-').find { it =~ /out .*\.(fa|fsa|fasta)/ } ?: 'out.fasta').replace('out ', '').trim()
    outinner_name   = (args.split('-').find { it =~ /outinner .*\.(fa|fsa|fasta)/ } ?: 'outinner.fasta').replace('outinner ', '').trim()

    def touch_out   = args.contains("-out")         ? "touch $out_name"             : ''
    def touch_inner = args.contains("-outinner")    ? "touch $outinner_name"        : ''
    """
    touch "${prefix}.${extension}"
    $touch_out
    $touch_inner

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(gt --version | head -1 | sed 's/gt (GenomeTools) //')
    END_VERSIONS
    """
}
