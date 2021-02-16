// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK_CREATESEQUENCEDICTIONARY {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::gatk4=4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.1.9.0--py39_0"
    }

    input:
    path fasta

    output:
    path "*.dict"        , emit: dict 
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CreateSequenceDictionary \\
        --REFERENCE $fasta \\
        $options.args

    echo \$(gatk CreateSequenceDictionary --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
