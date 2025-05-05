process BAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hdcf5f25_2' :
        'biocontainers/bamtools:2.5.2--hdcf5f25_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.{bed,fasta,fastq,json,pileup,sam,yaml}"), emit: data
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    format_cmd     = args ==~ /-format (bed|fasta|fastq|json|pileup|sam|yaml)/
    matched_format = args =~  /-format ([a-z]+)/
    extension      = matched_format[0][1]
    if ( format_cmd == false ) error "-format option must be provided in args. Possible values: bed fasta fastq json pileup sam yaml"

    """
    bamtools \\
        convert \\
        $args \\
        -in $bam \\
        -out ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    format_cmd     = args ==~ /-format (bed|fasta|fastq|json|pileup|sam|yaml)/
    matched_format = args =~  /-format ([a-z]+)/
    extension      = matched_format[0][1]
    if ( format_cmd == false ) error "-format option must be provided in args. Possible values: bed fasta fastq json pileup sam yaml"

    """
    touch ${prefix}.${extension}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
