process BAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c1/c1fde4cd562bddffee839570552212cbc3e01bb253ae8af6954a069106d518ba/data'
        : 'community.wave.seqera.io/library/bamtools:2.5.2--b804fdcee08dfdd6' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.{bed,fasta,fastq,json,pileup,sam,yaml}"), emit: data
    tuple val("${task.process}"), val('bamtools'), eval("bamtools --version | sed '2!d;s/bamtools //g'"), emit: versions_bamtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    format_cmd     = args ==~ /-format (bed|fasta|fastq|json|pileup|sam|yaml)/
    matched_format = args =~  /-format ([a-z]+)/
    extension      = matched_format[0][1]
    if ( format_cmd == false ) error "-format option must be provided in args. Possible values: bed fasta fastq json pileup sam yaml"

    """
    bamtools \\
        convert \\
        $args \\
        -in $bam \\
        -out ${prefix}.${extension}
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    format_cmd     = args ==~ /-format (bed|fasta|fastq|json|pileup|sam|yaml)/
    matched_format = args =~  /-format ([a-z]+)/
    extension      = matched_format[0][1]
    if ( format_cmd == false ) error "-format option must be provided in args. Possible values: bed fasta fastq json pileup sam yaml"

    """
    touch ${prefix}.${extension}
    """
}
