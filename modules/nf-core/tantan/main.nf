process TANTAN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tantan:50--h43eeafb_0':
        'biocontainers/tantan:50--h43eeafb_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fasta.gz") , optional:true , emit: masked_fasta
    tuple val(meta), path("${prefix}.tsv")      , optional:true , emit: repeat_probs
    tuple val(meta), path("${prefix}.tsv")      , optional:true , emit: repeat_counts
    tuple val(meta), path("${prefix}.bed")      , optional:true , emit: bed
    tuple val(meta), path("${prefix}.fasta.gz") , optional:true , emit: tandem_repeats
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def format_pattern  = /(-f)+\s+(\d)/
    def format_matcher  = (args =~ format_pattern)
    def format_num      = format_matcher[0][2]
    def output_format   = format_num == 1 || format_num == 2 ? ( format_num == 3 ? "bed" : "tsv" ) : "fasta"
    prefix              = task.ext.prefix ?: "${meta.id}"
    """
    tantan \\
        ${fasta} \\
        ${args} > ${prefix}.${output_format}

    if [ -f ${prefix}.fasta ]; then
        gzip ${prefix}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version 2>&1 | sed 's/^.*tantan //')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def format_pattern  = /(-f)+\s+(\d)/
    def format_matcher  = (args =~ format_pattern)
    def format_num      = format_matcher[0][2]
    def output_format   = format_num == 1 || format_num == 2 ? ( format_num == 3 ? "bed" : "tsv" ) : "fasta"
    prefix              = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${output_format}

    if [ -f ${prefix}.fasta ]; then
        gzip ${prefix}.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version 2>&1 | sed 's/^.*tantan //')
    END_VERSIONS
    """
}
