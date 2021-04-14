// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2_RUN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.1 conda-forge::pigz=2.6' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0'
    }

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path('*classified*')  , emit: classified
    tuple val(meta), path('*unclassified*'), emit: unclassified
    tuple val(meta), path('*report.txt')   , emit: txt
    path '*.version.txt'                   , emit: version

    script:
    def software     = getSoftwareName(task.process)
    def prefix       = options.suffix  ? "${meta.id}${options.suffix}"  : "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --unclassified-out $unclassified \\
        --classified-out $classified \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $paired \\
        $options.args \\
        $reads

    pigz -p $task.cpus *.fastq

    echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//' > ${software}.version.txt
    """
}
