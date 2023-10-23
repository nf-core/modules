process BLAST_BLASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0':
        'biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    val out_ext

    output:
    tuple val(meta), path("*.xml"), optional: true, emit: xml
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val(meta), path("*.csv"), optional: true, emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    switch ( out_ext ) {
        case "xml": outfmt = 5; break
        case "tsv": outfmt = 6; break
        case "csv": outfmt = 10; break
        default:
            outfmt = '6';
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting BLAST default of tabular BLAST output (tsv)");
            break
    }

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    DB=`find -L ./ -name "*.phr" | sed 's/\\.phr\$//'`
    blastp \\
        -query ${fasta_name} \\
        -out ${prefix}.${out_ext} \\
        -db \$DB \\
        -num_threads ${task.cpus} \\
        -outfmt ${outfmt} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    switch ( out_ext ) {
        case "xml": outfmt = 5; break
        case "tsv": outfmt = 6; break
        case "csv": outfmt = 10; break
        default:
            outfmt = '6';
            out_ext = 'tsv';
            log.warn("Unknown output file format provided (${out_ext}): selecting BLAST default of tabular BLAST output (tsv)");
            break
    }

    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """
}
