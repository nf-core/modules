process BLAST_BLASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
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
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
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
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.phr" | sed 's/\\.phr\$//'`
    blastp \\
        -query ${fasta_name} \\
        -out ${prefix}.${out_ext} \\
        -db \$DB \\
        -num_threads ${task.cpus} \\
        -outfmt ${outfmt} \\
        ${args}

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
