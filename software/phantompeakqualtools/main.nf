// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '1.2.2'

process PHANTOMPEAKQUALTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::phantompeakqualtools=1.2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/phantompeakqualtools:1.2.2--0"
    } else {
        container "quay.io/biocontainers/phantompeakqualtools:1.2.2--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.out")  , emit: spp
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.Rdata"), emit: rdata
    path  "*.version.txt"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    RUN_SPP=`which run_spp.R`
    Rscript -e "library(caTools); source(\\"\$RUN_SPP\\")" -c="$bam" -savp="${prefix}.spp.pdf" -savd="${prefix}.spp.Rdata" -out="${prefix}.spp.out" -p=$task.cpus
    echo $VERSION > ${software}.version.txt
    """
}
