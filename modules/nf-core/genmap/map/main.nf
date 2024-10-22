process GENMAP_MAP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1' :
        'biocontainers/genmap:1.3.0--h1b792b2_1' }"

    input:
    tuple val(meta), path(index)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.wig")      , optional:true, emit: wig
    tuple val(meta), path("*.bedgraph") , optional:true, emit: bedgraph
    tuple val(meta), path("*.txt")      , optional:true, emit: txt
    tuple val(meta), path("*.csv")      , optional:true, emit: csv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    def bed = regions ? "--selection ${regions}" : ""

    if ("$index" == "${prefix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    genmap \\
        map \\
        ${args} \\
        ${bed} \\
        --threads ${task.cpus} \\
        --index ${index} \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "$meta.id"
    def token_args = args.tokenize(" ")
    def wig = token_args.contains("-w") || token_args.contains("--wig") ?       "touch ${prefix}.wig"       : ""
    def bg =  token_args.contains("-bg") || token_args.contains("--bedgraph") ? "touch ${prefix}.bedgraph"  : ""
    def txt = token_args.contains("-t") || token_args.contains("--txt") ?       "touch ${prefix}.txt"       : ""
    def csv = token_args.contains("-d") || token_args.contains("--csv") ?       "touch ${prefix}.csv"       : ""

    if ("$index" == "${prefix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    ${wig}
    ${bg}
    ${txt}
    ${csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
