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
    tuple val(meta), path(tumourbam), path(normalbam)
    path fasta


    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.target.bed"), emit: bed
    tuple val(meta), path("*.antitarget.bed"), emit: bed
    tuple val(meta), path("*.targetcoverage.cnn"), emit: cnn
    tuple val(meta), path("*.antitargetcoverage.cnn"), emit: cnn
    tuple val(meta), path("*.cnn"), emit: cnn
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*.cns"), emit: cns
    tuple val(meta), path("*-scatter.pdf"), emit: pdf
    tuple val(meta), path("*-diagram.pdf"), emit: pdf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    if (meta.with_normal) {
	"""
        cnvkit.py batch $options.args $tumourbam \\ 
          --normal $normalbam \\
<<<<<<< HEAD
          --fasta $fasta \\
        cnvkit.py version > ${software}.version.txt
=======
          --method wgs \\
          --fasta $reffasta \\
          --annotate $annotationfile \\
          --output-reference reference.cnn --output-dir output
          cnvkit.py version > ${software}.version.txt
>>>>>>> e60433e152ce8d988c04b6a207debea4afb4326d
        """
    } else {
	"""
        cnvkit.py batch $options.args $tumourbam
          --normal \\
<<<<<<< HEAD
          --fasta $fasta \\
        cnvkit.py version > ${software}.version.txt
=======
          --method wgs \\
          --fasta $reffasta \\
          --annotate $annotationfile \\
          --output-reference my_flat_reference.cnn --output-dir output
          cnvkit.py version > ${software}.version.txt
>>>>>>> e60433e152ce8d988c04b6a207debea4afb4326d
        """
    }
}
