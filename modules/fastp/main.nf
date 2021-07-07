// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::fastp=0.20.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    } else {
        container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.json')         , emit: json
    tuple val(meta), path('*.html')         , emit: html
    tuple val(meta), path('*.log')          , emit: log
    path '*.version.txt'                    , emit: version
    tuple val(meta), path('*.fail.fastq.gz'), optional:true, emit: reads_fail

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        def fail_fastq = params.save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1 ${prefix}.trim.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $options.args \\
            2> ${prefix}.fastp.log
        echo \$(fastp --version 2>&1) | sed -e "s/fastp //g" > ${software}.version.txt
        """
    } else {
        def fail_fastq = params.save_trimmed_fail ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.trim.fastq.gz \\
            --out2 ${prefix}_2.trim.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $options.args \\
            2> ${prefix}.fastp.log

        echo \$(fastp --version 2>&1) | sed -e "s/fastp //g" > ${software}.version.txt
        """
    }
}
