// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def options    = initOptions(params.options)

process BEDTOOLS_MERGE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bedtools =2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"
    } else {
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    }
    input:
        path(sort)

    output:
        path("*.merged.bed"), emit: merge
        path  "*.version.txt", emit: version
// TODO fix output file naming
    script:
        def software = getSoftwareName(task.process)
        """
        bedtools merge -i $sort ${options.args} > test.merged.bed 
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
}
