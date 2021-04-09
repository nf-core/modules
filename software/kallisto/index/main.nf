// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    path fasta

    output:
    path "kallisto" , emit: idx
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    kallisto \\
        index \\
        $options.args \\
        -i kallisto \\
        $fasta

    echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > ${software}.version.txt
    """
}
