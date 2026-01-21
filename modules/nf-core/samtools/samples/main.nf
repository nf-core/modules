process SAMTOOLS_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0':
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta) , path(bam)  , path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Wrapping fasta in a list in case there is exactly one (in that case it's a bare path)
    def fasta_arg = fasta ? [fasta].flatten().collect { "-f $it" }.join(' ') : ''
    def out_arg = "-o ${prefix}.tsv"
    def bai_arg = args.contains('-X') ? bai : ''
    """
    samtools samples \\
        $args \\
        $fasta_arg \\
        $out_arg \\
        $bam \\
        $bai_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Write header if requested because the output will likely be parsed by Nextflow
    def headers = ''
    if (args.contains('-h\\n')) {
        headers = "#SM\tPATH"
        if (args.contains('-i')) {
            headers += "\tINDEX"
        }
        if (fasta) {
            headers += "\tREFERENCE"
        }
    }
    """
    echo $args

    echo -ne "$headers" > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
