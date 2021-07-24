// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path  '*.version.txt'         , emit: version
    tuple val(meta), path('*fastq.gz'), optional:true, emit: fastq

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        def unaligned = params.save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
        bowtie2 \\
            -x \$INDEX \\
            -U $reads \\
            --threads ${split_cpus} \\
            $unaligned \\
            $options.args \\
            2> ${prefix}.bowtie2.log \\
            | samtools view -@ ${split_cpus} $options.args2 -bhS -o ${prefix}.bam -

        echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
        """
    } else {
        def unaligned = params.save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
        bowtie2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --threads ${split_cpus} \\
            $unaligned \\
            $options.args \\
            2> ${prefix}.bowtie2.log \\
            | samtools view -@ ${split_cpus} $options.args2 -bhS -o ${prefix}.bam -

        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
        fi
        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
        fi
        echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
        """
    }
}
