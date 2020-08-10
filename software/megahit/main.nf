// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process MEGAHIT {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"

    conda (params.conda ? "bioconda::megahit=1.2.9-0" : null)

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    val(options)

    output:
    tuple val(meta), path("*.contigs.fa"), emit: fasta
    tuple val(meta), path("log"), emit: log

    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def inputs = ""
    if (meta.single_end) {
        inputs = "-r ${reads_1.join(',')}"
    } else if (meta.interleaved) {
        inputs = "--12 ${reads_1.join(',')}"
    } else {
        inputs = "-1 ${reads_1.join(',')} -2 ${reads_2.join(',')}"
    }
    """
    megahit $ioptions.args --num-cpu-threads $task.cpus $inputs -o results
    mv results/* .
    echo \$(megahit --version 2>&1) | sed 's/^.*MEGAHIT //; s/Using.*\$//' > ${software}.version.txt
    """
}
