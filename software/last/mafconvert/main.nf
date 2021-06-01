// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_MAFCONVERT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::last=1219" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1219--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1219--h2e03b76_0"
    }

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.{axt,blast,blasttab,chain,gff,html,psl,sam,tab}.gz"), emit: alm
    path "*.version.txt"                                                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def format   = params.options.format ? params.options.format  : "tab"
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    zcat $maf | maf-convert $options.args $format | gzip --no-name \\
        > ${prefix}.${format}.gz  

    # maf-convert has no --version option but lastdb (part of the same package) has.
    echo \$(lastdb --version 2>&1) | sed 's/lastdb //' > ${software}.version.txt
    """
}
