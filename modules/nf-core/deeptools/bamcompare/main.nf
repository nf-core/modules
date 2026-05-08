// This module was written employing Seqera AI (https://seqera.io/ask-ai/chat-v2)
process DEEPTOOLS_BAMCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0' :
        'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam1), path(bai1), path(bam2), path(bai2)

    output:
    tuple val(meta), path("*.bigWig")  , emit: bigwig  , optional: true
    tuple val(meta), path("*.bedgraph"), emit: bedgraph, optional: true
    tuple val("${task.process}"), val('deeptools'), eval('bamCompare --version | sed "s/bamCompare //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Determine output format from args or default to bigwig
    def format = args.contains('--outFileFormat bedgraph') ? 'bedgraph' : 'bigWig'
    def output_file = "${prefix}.${format}"

    """
    bamCompare \\
        --bamfile1 $bam1 \\
        --bamfile2 $bam2 \\
        --outFileName $output_file \\
        --numberOfProcessors $task.cpus \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = args.contains('--outFileFormat bedgraph') ? 'bedgraph' : 'bigWig'
    def output_file = "${prefix}.${format}"

    """
    touch $output_file

    """
}
