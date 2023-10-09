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
    tuple val(meta3), path(vcf)
    path fasta
    path fai

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
    calls_cmd = args.contains("-call") ? ">${prefix} + "_cnvnator.tab" : ''
    """
    cnvnator \\
        -root ${prefix}.pytor \\
        $args \\
        $reference \\
        $input_cmd \\
    	$calls_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
