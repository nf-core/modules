process BLAST_BLASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/4845ad74dd459a24f390f9929801119fcf93df10161cc65292e8f416ffc05682/data'
        : 'community.wave.seqera.io/library/blast:2.17.0--9196f264d54655d4' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val out_ext

    output:
    tuple val(meta), path("*.xml*") , emit: xml, optional: true
    tuple val(meta), path("*.tsv*") , emit: tsv, optional: true
    tuple val(meta), path("*.csv*") , emit: csv, optional: true
    tuple val("${task.process}"), val("blastp"), eval("blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//'"), topic: versions, emit: versions_blastp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''

    def outfmt = 6
    if ( "$out_ext" ==~ /^xml(\.gz)?$/ ) {
        outfmt = 5
    } else if ( "$out_ext" ==~ /^tsv(\.gz)?$/ ) {
        outfmt = 6
    } else if ( "$out_ext" ==~ /^csv(\.gz)?$/ ) {
        outfmt = 10
    } else {
        out_ext = 'tsv'
        outfmt = 6
        log.warn("Unknown output file format provided (${out_ext}): selecting BLAST default of tabular BLAST output (tsv)")
    }

    def out_ext_sans_gz = "$out_ext" - ~/\.gz$/
    def compress_output = "$out_ext" ==~ /^\w+\.gz$/ ? "gzip ${prefix}.${out_ext_sans_gz}" : ''

    """
    $uncompress_input

    DB=`find -L ./ -name "*.phr" | sed 's/\\.phr\$//'`
    blastp \\
        -query ${fasta_name} \\
        -out ${prefix}.${out_ext_sans_gz} \\
        -db \$DB \\
        -num_threads ${task.cpus} \\
        -outfmt ${outfmt} \\
        ${args}

    $compress_output

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def outfmt = 6
    if ( "$out_ext" ==~ /^xml(\.gz)?$/ ) {
        outfmt = 5
    } else if ( "$out_ext" ==~ /^tsv(\.gz)?$/ ) {
        outfmt = 6
    } else if ( "$out_ext" ==~ /^csv(\.gz)?$/ ) {
        outfmt = 10
    } else {
        out_ext = 'tsv'
        outfmt = 6
        log.warn("Unknown output file format provided (${out_ext}): selecting BLAST default of tabular BLAST output (tsv)")
    }

    def out_ext_sans_gz = "$out_ext" - ~/\.gz$/
    def compress_output = "$out_ext" ==~ /^\w+\.gz$/ ? "gzip ${prefix}.${out_ext_sans_gz}" : ''

    """
    touch ${prefix}.${out_ext_sans_gz}

    $compress_output

    """
}
