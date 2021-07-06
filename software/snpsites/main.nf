// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPSITES {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::snp-sites=2.5.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0"
    } else {
        container "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
    }

    input:
    path alignment

    output:
    path "*.fas"        , emit: fasta
    path "*.sites.txt"  , emit: constant_sites
    path "*.version.txt", emit: version
    env   CONSTANT_SITES, emit: constant_sites_string

    script:
    def software = getSoftwareName(task.process)
    """
    snp-sites \\
        $alignment \\
        $options.args \\
        > filtered_alignment.fas

    echo \$(snp-sites -C $alignment) > constant.sites.txt

    CONSTANT_SITES=\$(cat constant.sites.txt)

    echo \$(snp-sites -V 2>&1) | sed 's/snp-sites //' > ${software}.version.txt
    """
}
