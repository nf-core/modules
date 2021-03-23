// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fgbio:1.3.0--0"
    } else {
        container "quay.io/biocontainers/fgbio:1.3.0--0"
    }

    input:
    tuple val(meta), path(reads)
    val(read_structure1)
    val(read_structure2)

    output:
    tuple val(meta), path("*_umi_converted.bam"), emit: bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir tmpFolder

    fgbio --tmp-dir=${PWD}/tmpFolder \
    FastqToBam \
    -i ${reads} \
    -o "${prefix}_umi_converted.bam" \
    --read-structures $read_structure1 $read_structure2 \
    --sample ${meta.id} \
    --library ${meta.id}

    echo \$(fgbio --version 2>&1) >${software}.version.txt
    """
}
