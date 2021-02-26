// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CNVKIT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::cnvkit=0.9.8=0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/conda:0.9.8--py_0"
    } else {
        container "quay.io/biocontainers/cnvkit:0.9.8--py_0"
    }


    input:
    tuple val(meta), path(bam)
    path fasta


    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*.cns"), emit: cns
    tuple val(meta), path("*.pdf"), emit: pdf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_tumour_278_sub_chr21.bam ] && ln -s ${bam[0]} ${prefix}_tumour_278_sub_chr21.bam
    [ ! -f  ${prefix}_normal_280_sub_chr21.bam ] && ln -s ${bam[1]} ${prefix}_normal_280_sub_chr21.bam
    cnvkit.py batch $options.args ${prefix}_tumour_278_sub_chr21.bam\\
        --normal ${prefix}_normal_280_sub_chr21.bam\\
        --fasta $fasta \\
    cnvkit.py version > ${software}.version.txt
    """
}
