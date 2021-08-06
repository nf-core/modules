// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COOLER_DIGEST {
    tag '$fasta'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0"
    }

    input:
    path fasta
    path chromsizes
    val  enzyme

    output:
    path "*.bed"                  , emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    cooler digest \\
        $options.args \\
        -o "${fasta.baseName}_${enzyme.replaceAll(/[^0-9a-zA-Z]+/, "_")}.bed" \\
        $chromsizes \\
        $fasta \\
        $enzyme

    echo \$(cooler --version 2>&1) | sed 's/cooler, version //' > ${software}.version.txt
    """
}
