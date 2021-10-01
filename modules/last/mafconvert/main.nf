// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_MAFCONVERT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1250--h2e03b76_0"
    }

    input:
    tuple val(meta), path(maf)
    val(format)

    output:
    tuple val(meta), path("*.axt.gz"),      optional:true, emit: axt_gz
    tuple val(meta), path("*.blast.gz"),    optional:true, emit: blast_gz
    tuple val(meta), path("*.blasttab.gz"), optional:true, emit: blasttab_gz
    tuple val(meta), path("*.chain.gz"),    optional:true, emit: chain_gz
    tuple val(meta), path("*.gff.gz"),      optional:true, emit: gff_gz
    tuple val(meta), path("*.html.gz"),     optional:true, emit: html_gz
    tuple val(meta), path("*.psl.gz"),      optional:true, emit: psl_gz
    tuple val(meta), path("*.sam.gz"),      optional:true, emit: sam_gz
    tuple val(meta), path("*.tab.gz"),      optional:true, emit: tab_gz
    path "versions.yml"                                  , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    maf-convert $options.args $format $maf | gzip --no-name \\
        > ${prefix}.${format}.gz

    # maf-convert has no --version option but lastdb (part of the same package) has.
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
