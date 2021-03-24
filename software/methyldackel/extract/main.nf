// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::methyldackel=0.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/methyldackel:0.5.2--h7435645_0"
    } else {
        container "quay.io/biocontainers/methyldackel:0.5.2--h7435645_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path  "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    MethylDackel extract \\
        $options.args \\
        $fasta \\
        $bam

    echo \$(methyldackel --version 2>&1) | cut -f1 -d" " > ${software}.version.txt
    """
}
