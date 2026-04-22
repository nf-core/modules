process SAMTOOLS_SAMPLES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Wrapping fasta in a list in case there is exactly one (in that case it's a bare path)
    def fasta_arg = fasta ? [fasta].flatten().collect { fasta_file -> "-f ${fasta_file}" }.join(' ') : ''
    def out_arg = "-o ${prefix}.tsv"
    def bai_arg = args.contains('-X') ? bai : ''
    """
    samtools samples \\
        ${args} \\
        ${fasta_arg} \\
        ${out_arg} \\
        ${bam} \\
        ${bai_arg}
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
    echo ${args}
    echo -ne "${headers}" > ${prefix}.tsv
    """
}
