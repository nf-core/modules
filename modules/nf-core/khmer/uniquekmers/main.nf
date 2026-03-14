process KHMER_UNIQUEKMERS {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    tuple val(meta), path(fasta)
    val kmer_size

    output:
    tuple val(meta), path("*.report.txt")  , emit: report
    tuple val(meta), path("*.kmers.txt")   , emit: kmers
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    unique-kmers.py \\
        -k $kmer_size \\
        -R ${prefix}.report.txt \\
        $args \\
        $fasta

    grep ^number ${prefix}.report.txt | sed 's/^.*:.[[:blank:]]//g' > ${prefix}.kmers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( unique-kmers.py --version 2>&1 | grep ^khmer | sed 's/^khmer //;s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.txt
    touch ${prefix}.kmers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( unique-kmers.py --version 2>&1 | grep ^khmer | sed 's/^khmer //;s/ .*\$//' )
    END_VERSIONS
    """
}
