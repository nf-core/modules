process ULTRAPLEX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ultraplex=1.2.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ultraplex:1.2.5--py38h4a8c8d9_0' :
        'biocontainers/ultraplex:1.2.5--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(reads)
    path(barcode_file)
    val(adapter_seq)

    output:
    tuple val(meta), path("*[!no_match].fastq.gz"),              emit: fastq
    tuple val(meta), path("*no_match.fastq.gz"), optional: true, emit: no_match_fastq
    path "*.log",                                                emit: report
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "1.2.5" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    def adapter_seq_command = ''
    if(adapter_seq) {
        adapter_seq_command = "--adapter ${adapter_seq}"
    }

    read_list = reads.collect{it.toString()}
    if (read_list.size > 1){
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --input_2 ${read_list[1]} \\
        --barcodes $barcode_file \\
        --threads $task.cpus $args $adapter_seq_command"""
    } else {
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --barcodes $barcode_file \\
        --threads $task.cpus $args $adapter_seq_command"""
    }

    """
    ${ultraplex_command}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultraplex: $VERSION
    END_VERSIONS
    """
}
