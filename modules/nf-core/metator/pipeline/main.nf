process METATOR_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46a2ce37dac94fcc2920f7b4d42c60455fe26b00945ff5e51624cba61373e8f9/data'
        : 'community.wave.seqera.io/library/metator_pysam_python:eaa4c24c77402ae1'}"

    input:
    tuple val(meta), path(contigs), path(hic_input, arity: '1..2'), path(depths)

    output:
    tuple val(meta), path("final_bin_unscaffold/*.fa.gz"), emit: bins
    tuple val(meta), path("bin_summary.txt")             , emit: bin_summary
    tuple val(meta), path("binning.txt")                 , emit: contig2bin
    tuple val(meta), path("network.txt")                 , emit: network
    tuple val(meta), path("contig_data_final.txt")       , emit: contig_data
    tuple val(meta), path("plot/*.png")                  , emit: plots, optional: true
    tuple val("${task.process}"), val('metator'), eval("metator -v"), topic: versions, emit: versions_metator
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version |& sed "1!d;s/^.*(gzip) //;s/ Copyright.*//"'), topic: versions, emit: versions_gunzip
    tuple val("${task.process}"), val("find"), eval("find --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_find

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def depth_input = depths ? "--depth ${depths}" : ""
    def assembly_input = contigs =~ /\.gz$/ ? "${contigs.getBaseName()}" : contigs
    def gunzip = contigs =~ /\.gz$/ ? "gunzip ${args2} -c ${contigs} > ${assembly_input}" : ""

    // Set up type input
    def input_type_arg = "-S fastq"
    if (hic_input[0].getExtension() == "bam") {
        input_type_arg = "-S bam"
    } else if (hic_input[0].getName().endsWith(".pairs.gz") || hic_input[0].getName().endsWith(".pairs")) {
        input_type_arg = "-S pair"
        if (hic_input.size() > 1) {
            error("Error: Pairs file supplied to Metator but a second Hi-C file also supplied. Please only supply a pairs file.")
        }
    }

    // Set up input file args
    def input_arg = "--forward ${hic_input[0]}" +
        (hic_input.size() == 2 ? " --reverse ${hic_input[1]}" : "")
    """
    ${gunzip}

    metator pipeline \\
        ${input_arg} \\
        --assembly ${assembly_input} \\
        ${depth_input} \\
        ${input_type_arg} \\
        -t ${task.cpus} \\
        --prefix ${prefix} \\
        ${args}

    find \\
        final_bin_unscaffold/ \\
        ${args3} \\
        -name "*.fa" \\
        -exec gzip {} \\;
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch bin_summary.txt
    touch binning.txt
    touch network.txt
    touch contig_data_final.txt

    mkdir plot
    touch plot/bins_distribution.png
    touch plot/bins_size_distribution.png
    touch plot/MAGs_GC_distribution.png
    touch plot/MAGs-HiC_cov_distribution.png

    mkdir final_bin_unscaffold
    echo "" | gzip > final_bin_unscaffold/${prefix}_metator_00001_00000.fa.gz
    """
}
