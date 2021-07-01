// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_POSTMASK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::last=1238" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1238--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1238--h2e03b76_0"
    }

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use the suffix option to disambiguate"
    """
    last-postmask $options.args $maf | gzip --no-name > ${prefix}.maf.gz

    # last-postmask does not have a --version option
    echo \$(lastal --version 2>&1) | sed 's/lastal //' > ${software}.version.txt
    """
}
