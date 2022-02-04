process BAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamtools=2.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.1--h9a82719_9' :
        'quay.io/biocontainers/bamtools:2.5.1--h9a82719_9' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.{bed,fasta,fastq,json,pileup,sam,yaml}"), emit: data
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def test = args ==~ /-format (bed|fasta|fastq|json|pileup|sam|yaml)/
    if ( test == false ) error "-format option must be provided in args. Possible values: bed fasta fastq json pileup sam yaml"
    m = args =~ /-format ([a-z]+)/
    ext = m[0][1]


    """
    bamtools \\
        convert \\
        $args \\
        -in $bam \\
        -out ${prefix}.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
