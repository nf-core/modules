process DIAMOND_BLASTP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.11--h5ca1c30_0' :
        'biocontainers/diamond:2.1.11--h5ca1c30_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val outfmt
    val blast_columns

    output:
    tuple val(meta), path('*.{blast,blast.gz}'), optional: true, emit: blast
    tuple val(meta), path('*.{xml,xml.gz}')    , optional: true, emit: xml
    tuple val(meta), path('*.{txt,txt.gz}')    , optional: true, emit: txt
    tuple val(meta), path('*.{daa,daa.gz}')    , optional: true, emit: daa
    tuple val(meta), path('*.{sam,sam.gz}')    , optional: true, emit: sam
    tuple val(meta), path('*.{tsv,tsv.gz}')    , optional: true, emit: tsv
    tuple val(meta), path('*.{paf,paf.gz}')    , optional: true, emit: paf
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def columns = blast_columns ? "${blast_columns}" : ''
    def out_ext = ""

    if (outfmt == 0) {
        out_ext = "blast"
    } else if (outfmt == 5) {
        out_ext = "xml"
    } else if (outfmt == 6) {
        out_ext = "txt"
    } else if (outfmt == 100) {
        out_ext = "daa"
    } else if (outfmt == 101) {
        out_ext = "sam"
    } else if (outfmt == 102) {
        out_ext = "tsv"
    } else if (outfmt == 103) {
        out_ext = "paf"
    } else {
        outfmt = '6'
        out_ext = 'txt'
        log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)")
    }

    if ( args =~ /--compress\s+1/ ) out_ext += '.gz'

    """
    diamond \\
        blastp \\
        --threads ${task.cpus} \\
        --db ${db} \\
        --query ${fasta} \\
        --outfmt ${outfmt} ${columns} \\
        ${args} \\
        --out ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def out_ext = ""

    if (outfmt == 0) {
        out_ext = "blast"
    } else if (outfmt == 5) {
        out_ext = "xml"
    } else if (outfmt == 6) {
        out_ext = "txt"
    } else if (outfmt == 100) {
        out_ext = "daa"
    } else if (outfmt == 101) {
        out_ext = "sam"
    } else if (outfmt == 102) {
        out_ext = "tsv"
    } else if (outfmt == 103) {
        out_ext = "paf"
    } else {
        outfmt = '6'
        out_ext = 'txt'
        log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)")
    }

    if ( args =~ /--compress\s+1/ ) out_ext += '.gz'

    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
