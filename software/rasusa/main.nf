// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RASUSA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::rasusa=0.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rasusa:0.3.0--h779adbc_1"
    } else {
        container "quay.io/biocontainers/rasusa:0.3.0--h779adbc_1"
    }

    input:
    tuple val(meta), path(reads)
    val depth_cutoff
    val genome_size

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path '*.version.txt'                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}_${depth_cutoff}X" : "${meta.id}_${depth_cutoff}X"
    def output = meta.single_end ? "--output ${prefix}.fastq.gz" : "--output ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz" 
    """
    rasusa \\
        $options.args \\
        --coverage $depth_cutoff \\
        --genome-size $genome_size \\
        --input $reads \\
        $output
    echo \$(rasusa --version 2>&1) | sed -e "s/rasusa //g" > ${software}.version.txt
    """
}
