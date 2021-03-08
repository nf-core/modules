// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::bcftools=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.11--h7c999a4_0'
    }

    input:
    tuple val(meta), path(vcf), path(tbi), path(fasta)

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    path  '*.version.txt'        , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cat $fasta | bcftools consensus $vcf $options.args > ${prefix}.fa
    header=\$(head -n 1 ${prefix}.fa | sed 's/>//g')
    sed -i 's/\${header}/${meta.id}/g' ${prefix}.fa

    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
