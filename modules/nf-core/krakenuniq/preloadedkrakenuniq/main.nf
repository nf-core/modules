process KRAKENUNIQ_PRELOADEDKRAKENUNIQ {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::krakenuniq=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.0--pl5321h19e8d03_0':
        'quay.io/biocontainers/krakenuniq:1.0.0--pl5321h19e8d03_0' }"

    input:
    tuple val(meta), path(fastqs)
    path  db
    val ram_chunk_size
    val save_output_fastqs
    val report_file
    val save_output

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classified.txt')        , optional:true, emit: classified_assignment
    tuple val(meta), path('*report.txt')                           , emit: report

    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''

    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? '"\$PREFIX".classified.fastq'   : '"\$PREFIX".classified#.fastq'
    def unclassified = meta.single_end ? '"\$PREFIX".unclassified.fastq' : '"\$PREFIX".unclassified#.fastq'
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def output_option = save_output ? '--output "\$PREFIX".krakenuniq.classified.txt' : ""
    def report = report_file ? '--report-file "\$PREFIX".krakenuniq.report.txt' : ""
    def compress_reads_command = save_output_fastqs ? "gzip *.fastq" : ""

    """
    krakenuniq \\
        $args \\
        --db $db \\
        --preload $ram_chunk_size \\
        --threads $task.cpus;

    for fastq in ${fastqs.join(' ')}; do \\
        PREFIX=\$(echo \$fastq);
        krakenuniq \\
            --db $db \\
            --threads $task.cpus \\
            $report \\
            $output_option \\
            $unclassified_option \\
            $classified_option \\
            $output_option \\
            $paired \\
            $args2 \\
            \$fastq;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    echo stub > test_1.fastq.gz.classified.fastq
    echo stub > test_1.fastq.gz.krakenuniq.classified.txt
    echo stub > test_1.fastq.gz.krakenuniq.report.txt
    echo stub > test_1.fastq.gz.unclassified.fastq
    echo stub > test_2.fastq.gz.classified.fastq
    echo stub > test_2.fastq.gz.krakenuniq.classified.txt
    echo stub > test_2.fastq.gz.krakenuniq.report.txt
    echo stub > test_2.fastq.gz.unclassified.fastq

    echo "${task.process}:" > versions.yml
    echo ' krakenuniq: 1.0.0' >> versions.yml
    """
}
