// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process CIRIQUANT_QUANT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.3--pyhdfd78af_0':
        'quay.io/biocontainers/ciriquant' }"

    input:
    tuple val(meta), path(reads)
    path(config)
    path(bed) //optional
    tuple val(meta1), path(bam) //optional
    tuple val(meta2), path(circ) //optional
    tuple val(meta3), path(rnaser) //optional

    output:
    //
    tuple val(meta), path("*.gtf"), emit: gtf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("CIRIerror.log"), emit: error_log, optional: true

    // Reference alignment outputs

    tuple val(meta), path("align/*.sorted.bam"), emit: sorted_bam
    tuple val(meta), path("align/*.sorted.bam.bai"), emit: sorted_bai

    // Circ RNA detection

    tuple val(meta), path("circ/*.ciri"), emit: ciri
    tuple val(meta), path("circ/*.ciri.bed"), emit: ciri_bed

    // De Novo Back Splice Junction Alignment Tracks

    tuple val(meta), path("circ/*_denovo.sorted.bam"), emit: denovo_sorted_bam
    tuple val(meta), path("circ/*_denovo.sorted.bam.bai"), emit: denovo_sorted_bai

    // Version Broadcaste

    tuple val("${task.process}"), val('ciriquant'), eval("CIRIquant --version"), topic: versions, emit: versions_ciriquant

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Optional Inputs

    def bed_option = bed ? "--bed ${bed}" : ""
    def bam_option = bam ? "--bam ${bam[1]}" : ""
    def circ_option = circ ? "--circ ${circ[1]}" : ""
    def rnaser_option = rnaser ? "--RNaseR ${rnaser[1]}" : ""

    """
    ciriquant \\
        $args \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --config $config \\
        -o ./ \\
        -p $prefix \\
        $bed_option \\
        $bam_option \\
        $circ_option \\
        $rnaser_option
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Create sub directories expected by ciriquant
    """
    mkdir -p align circ

    touch ${prefix}.gtf
    touch ${prefix}.log
    touch ${prefix}.bed
    touch ${prefix}.CIRIerror.log

    touch align/${prefix}.sorted.bam
    touch align/${prefix}.sorted.bam.bai

    touch circ/${prefix}.ciri
    touch circ/${prefix}.ciri.bed

    touch circ/${prefix}_denovo.sorted.bam
    touch circ/${prefix}__denovo.sorted.bam.bai
    """
}
