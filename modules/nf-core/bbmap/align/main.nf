// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(fastq)
    path ref
    tuple val(output_format), val(output_reads), val(interleaved_output)

    output:
    tuple val(meta), path("*.bam")                            , optional: true, emit: bam
    tuple val(meta), path("*.{fastq,fastq.gz,fasta,fasta.gz}"), optional: true, emit: reads
    tuple val(meta), path("*.log")                            , emit: log
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = bbmapAlignOutputs(meta, prefix, output_format, output_reads, interleaved_output)

    def fastq_files = fastq instanceof List ? fastq : [fastq]
    def input = fastq_files.size() == 2 ? "in=${fastq_files[0]} in2=${fastq_files[1]}" :
        meta.single_end ? "in=${fastq_files[0]}" : "in=${fastq_files[0]} interleaved=t"
    def output = outputs.files.size() == 2 ?
        "${outputs.key}=${outputs.files[0]} ${outputs.key}2=${outputs.files[1]}" :
        "${outputs.key}=${outputs.files[0]}"
    // An output named like a staged input makes bbmap fail, or clobbers the input with overwrite=t
    if ( outputs.files.any { file_name -> file_name in fastq_files*.name } ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }
    // Output channels are optional, so check bbmap wrote them; only uncompressed fastq/fasta may be empty
    def check_outputs = outputs.files.collect { file_name ->
        def test_flag = file_name.endsWith('.bam') || file_name.endsWith('.gz') ? '-s' : '-e'
        "test ${test_flag} ${file_name} || { echo 'ERROR: bbmap did not write ${file_name}' >&2; exit 1; }"
    }.join('\n    ')

    // Set the db variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containing a 'ref' subdir) or 3) a sequence
    // file in fasta format
    def db = ''
    if ( ref.isDirectory() ) {
        if ( !(ref ==~ /(.\/)?ref\/?/) ) {
            db = "path=${ref}"
        }
    } else {
        db = "ref=${ref}"
    }

    """
    bbmap.sh \\
        $db \\
        $input \\
        $output \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbmap.log

    ${check_outputs}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = bbmapAlignOutputs(meta, prefix, output_format, output_reads, interleaved_output)
    def create_outputs = outputs.files.collect { file_name ->
        file_name.endsWith('.gz') ? "echo '' | gzip > ${file_name}" : "touch ${file_name}"
    }.join('\n    ')
    """
    ${create_outputs}
    touch ${prefix}.bbmap.log
    """
}

// Shared by script and stub so both validate and name outputs identically
def bbmapAlignOutputs(meta, prefix, output_format, output_reads, interleaved_output) {
    def format = output_format ?: 'bam'
    def reads_to_write = output_reads ?: 'all'
    def valid_formats = ['bam', 'fastq', 'fastq.gz', 'fasta', 'fasta.gz']
    def valid_reads = ['all', 'mapped', 'unmapped']
    if ( !(format in valid_formats) ) {
        error "Invalid output_format '${format}'. Valid options: ${valid_formats.join(', ')}"
    }
    if ( !(reads_to_write in valid_reads) ) {
        error "Invalid output_reads '${reads_to_write}'. Valid options: ${valid_reads.join(', ')}"
    }
    // Unsuffixed only for all reads as bam, the module's original output
    def out_prefix = format == 'bam' && reads_to_write == 'all' ? prefix : "${prefix}.${reads_to_write}"
    def split_pairs = !meta.single_end && !interleaved_output && format != 'bam'
    def files = split_pairs ?
        [ "${out_prefix}_1.${format}", "${out_prefix}_2.${format}" ] :
        [ "${out_prefix}.${format}" ]
    return [
        key  : [ all: 'out', mapped: 'outm', unmapped: 'outu' ][reads_to_write],
        files: files.collect { file_name -> file_name.toString() }
    ]
}
