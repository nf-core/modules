process CNVNATOR_CNVNATOR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::cnvnator=0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvnator:0.4.1--py310h2dce045_7':
        'biocontainers/cnvnator:0.4.1--py310h2dce045_7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(pytor)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(output_meta), path("${prefix}.pytor")       , emit: pytor
    tuple val(output_meta), path("${prefix}_cnvnator.tab"), emit: tab, optional: true
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def input_cmd = bam             ? "-tree ${bam}"      : ''
    output_meta   = bam             ? meta                : meta2
    prefix        = task.ext.prefix ?: bam ? "${meta.id}" : "${meta2.id}"
    if (fasta) {
        reference = fasta.isDirectory() ? "-d ${fasta}" : "-fasta ${fasta}"
    } else {
        reference = ''
    }
    calls_cmd = args.contains("-call") ? ">${prefix}_cnvnator.tab" : ''
    """
    cnvnator \\
        -root ${prefix}.pytor \\
        $args \\
        $reference \\
        $input_cmd \\
    	$calls_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(cnvnator 2>&1) | sed -n '/CNVnator/p' | sed 's/CNVnator v//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def calls_cmd = args.contains("-call") ? "touch ${prefix}_cnvnator.tab" : ''
    """
    touch ${prefix}.pytor
    $calls_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(cnvnator 2>&1) | sed -n '/CNVnator/p' | sed 's/CNVnator v//')
    END_VERSIONS
    """
}
