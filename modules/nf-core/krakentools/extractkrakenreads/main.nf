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
    tuple val(meta), path("*.{fastq,fasta}"), emit: extracted_kraken2_reads

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads_command = meta.single_end ? "-s $classified_reads_fastq" : "-s1 ${classified_reads_fastq[0]} -s2 ${classified_reads_fastq[1]}"
    def fastq_out = task.ext.fastq_output ? "--fastq-output" : ""
    def output_format = task.ext.fastq_output ? "fastq" : "fasta"
    def output_reads_command = meta.single_end ? "-o ${prefix}.extracted_kraken2_read.$output_format" : "-o ${prefix}.extracted_kraken2_read_1.$output_format -o2 ${prefix}.extracted_kraken2_read_2.$output_format"
    def include_children_option = task.ext.include_children ? "--include-children" : ""
    def include_parents_option = task.ext.include_parents ? "--include-parents" : ""
    def report_option = (task.ext.include_children || task.ext.include_parents) ? "-r ${report}" : ""
    def exclude_option = task.ext.exlude ? "--exclude" : ""
    def append_option = task.ext.append ? "--append" : "--noappend"
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    extract_kraken_reads.py \\
        -t $taxid \\
        -k $classified_reads_assignment \\
        $input_reads_command \\
        $output_reads_command \\
        $report_option \\
        $include_children_option \\
        $include_parents_option \\
        $fastq_out \\
        $append_option \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """
}
