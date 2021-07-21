// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BISMARK_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    path fasta, stageAs: "BismarkIndex/*"

    output:
    path "BismarkIndex" , emit: index
    path "*.version.txt", emit: version

    script:
    def software   = getSoftwareName(task.process)
    """
    bismark_genome_preparation \\
        $options.args \\
        BismarkIndex

    echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//' > ${software}.version.txt
    """
}
