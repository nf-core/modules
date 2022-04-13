process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.2 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_classified
    val save_readclassification

    output:
    tuple val(meta), path('*classified*')     , optional:true, emit: classified
    tuple val(meta), path('*unclassified*')   , optional:true, emit: unclassified
    tuple val(meta), path('*classifiedreads*'), optional:true, emit: classifiedreads
    tuple val(meta), path('*report.txt')                     , emit: txt
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_command = save_classified ? "--classified-out ${classified}" : ""
    def unclassified_command = save_classified ? "--unclassified-out ${unclassified}" : ""
    def readclassification_command = save_readclassification ? "--output ${prefix}.kraken2.classifiedreads.txt" : ""
    def compress_reads_command = save_classified ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $unclassified_command \\
        $classified_command \\
        $readclassification_command \\
        $paired \\
        $args \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
