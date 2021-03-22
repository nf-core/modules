// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::cancerit-allelecount=4.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.2.1--h3ecb661_0"
    } else {
        container "quay.io/biocontainers/cancerit-allelecount:4.2.1--h3ecb661_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path loci

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    path "*.version.txt"                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    alleleCounter \\
        $options.args \\
        -l $loci \\
        -b $bam \\
        -o ${prefix}.alleleCount

    alleleCounter --version > ${software}.version.txt
    """
}
