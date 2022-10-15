process KHMER_UNIQUEKMERS {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    def container_image = "khmer:3.0.0a3--py37haa7609a_2"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    path fasta
    val  kmer_size

    output:
    path "report.txt"  , emit: report
    path "kmers.txt"   , emit: kmers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    unique-kmers.py \\
        -k $kmer_size \\
        -R report.txt \\
        $args \\
        $fasta

    grep ^number report.txt | sed 's/^.*:.[[:blank:]]//g' > kmers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( unique-kmers.py --version 2>&1 | grep ^khmer | sed 's/^khmer //;s/ .*\$//' )
    END_VERSIONS
    """
}
