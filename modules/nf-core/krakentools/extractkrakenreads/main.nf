process KRAKENTOOLS_EXTRACTKRAKENREADS {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0'}"

    input:
    val taxid // Separated by spaces
    tuple val(meta), path(classified_reads_assignment)
    tuple val(meta), path(classified_reads_fastq)
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.{fastq.gz,fasta.gz}"), emit: extracted_kraken2_reads

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--fastq-output") ? "fastq" : "fasta"
    def input_reads_command = meta.single_end ? "-s $classified_reads_fastq" : "-s1 ${classified_reads_fastq[0]} -s2 ${classified_reads_fastq[1]}"
    def output_reads_command = meta.single_end ? "-o ${prefix}.extracted_kraken2_read.${extension}" : "-o ${prefix}.extracted_kraken2_read_1.${extension} -o2 ${prefix}.extracted_kraken2_read_2.${extension}"
    def gzip_reads_command = meta.single_end ? "gzip ${prefix}.extracted_kraken2_read.${extension}" : "gzip ${prefix}.extracted_kraken2_read_1.${extension}; gzip ${prefix}.extracted_kraken2_read_2.${extension}"
    def report_option = report ? "-r ${report}" : ""
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    extract_kraken_reads.py \\
        ${args} \\
        -t $taxid \\
        -k $classified_reads_assignment \\
        $report_option \\
        $input_reads_command \\
        $output_reads_command

    $gzip_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--fastq-output") ? "fastq" : "fasta"
    def input_reads_command = meta.single_end ? "-s $classified_reads_fastq" : "-s1 ${classified_reads_fastq[0]} -s2 ${classified_reads_fastq[1]}"
    def output_reads_command = meta.single_end ? "-o ${prefix}.extracted_kraken2_read.${extension}" : "-o ${prefix}.extracted_kraken2_read_1.${extension} -o2 ${prefix}.extracted_kraken2_read_2.${extension}"
    def gzip_reads_command = meta.single_end ? "gzip ${prefix}.extracted_kraken2_read.${extension}" : "gzip ${prefix}.extracted_kraken2_read_1.${extension}; gzip ${prefix}.extracted_kraken2_read_2.${extension}"
    def report_option = report ? "-r ${report}" : ""
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    if [ "$meta.single_end" == "true" ];
    then
        touch ${prefix}.extracted_kraken2_read.${extension}
        $gzip_reads_command
    else
        touch ${prefix}.extracted_kraken2_read_1.${extension}
        touch ${prefix}.extracted_kraken2_read_2.${extension}
        $gzip_reads_command
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """
}
