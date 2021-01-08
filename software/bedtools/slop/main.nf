// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
def options    = initOptions(params.options)

process BEDTOOLS_SLOP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bedtools =2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"
    } else {
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    }


    input:
        tuple val(meta), path(beds), path (sizes)

    output:
        tuple val(meta), path("*.slop.bed"), emit: slopbed
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        def header = params.header ? "-header":''
        def pct = params.pct ? "-pct":''

        def symmetry = ''
        if (meta.symmetry) {

        """
        slopBed -i $beds -g $sizes -b $params.b $header $pct $options.args> ${prefix}.slop.bed
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
        } else {
        """
        slopBed -i $beds -g $sizes -l $params.l -r $params.r $header $pct $options.args> ${prefix}.slop.bed
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
        }
}
