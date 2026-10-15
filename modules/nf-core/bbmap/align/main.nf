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
    val output_format
    val output_reads
    val interleaved

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

    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    def out_key = [ all: 'out', mapped: 'outm', unmapped: 'outu' ][reads_to_write]
    def out_prefix = reads_to_write == 'all' ? prefix : "${prefix}.${reads_to_write}"
    def output_files = (meta.single_end || interleaved || format == 'bam') ?
        [ "${out_prefix}.${format}" ] :
        [ "${out_prefix}_1.${format}", "${out_prefix}_2.${format}" ]
    def output = output_files.size() == 1 ?
        "${out_key}=${output_files[0]}" :
        "${out_key}=${output_files[0]} ${out_key}2=${output_files[1]}"
    // An output named like a staged input makes bbmap fail, or clobbers the input with overwrite=t
    def input_names = "${fastq}".tokenize(' ')
    if ( output_files.any { file_name -> input_names.contains(file_name.toString()) } ) {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = output_format ?: 'bam'
    def reads_to_write = output_reads ?: 'all'
    def out_prefix = reads_to_write == 'all' ? prefix : "${prefix}.${reads_to_write}"
    def output_files = (meta.single_end || interleaved || format == 'bam') ?
        [ "${out_prefix}.${format}" ] :
        [ "${out_prefix}_1.${format}", "${out_prefix}_2.${format}" ]
    def create_outputs = output_files.collect { file_name ->
        file_name.endsWith('.gz') ? "echo '' | gzip > ${file_name}" : "touch ${file_name}"
    }.join('\n    ')
    """
    ${create_outputs}
    touch ${prefix}.bbmap.log
    """
}
