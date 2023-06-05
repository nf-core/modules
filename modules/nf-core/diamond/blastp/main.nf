process DIAMOND_BLASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::diamond=2.0.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_0' :
        'biocontainers/diamond:2.0.15--hb97b32f_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    val out_ext
    val blast_columns

    output:
    tuple val(meta), path('*.blast'), optional: true, emit: blast
    tuple val(meta), path('*.xml')  , optional: true, emit: xml
    tuple val(meta), path('*.txt')  , optional: true, emit: txt
    tuple val(meta), path('*.daa')  , optional: true, emit: daa
    tuple val(meta), path('*.sam')  , optional: true, emit: sam
    tuple val(meta), path('*.tsv')  , optional: true, emit: tsv
    tuple val(meta), path('*.paf')  , optional: true, emit: paf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = blast_columns ? "${blast_columns}" : ''
    switch ( out_ext ) {
        case "blast": outfmt = 0; break
        case "xml": outfmt = 5; break
        case "txt": outfmt = 6; break
        case "daa": outfmt = 100; break
        case "sam": outfmt = 101; break
        case "tsv": outfmt = 102; break
        case "paf": outfmt = 103; break
        default:
            outfmt = '6';
            out_ext = 'txt';
            log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)");
            break
    }
    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/\\.dmnd\$//'`

    diamond \\
        blastp \\
        --threads $task.cpus \\
        --db \$DB \\
        --query $fasta \\
        --outfmt ${outfmt} ${columns} \\
        $args \\
        --out ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
