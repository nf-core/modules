// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::methyldackel=0.6.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/methyldackel:0.6.0--h22771d5_0"
    } else {
        container "quay.io/biocontainers/methyldackel:0.6.0--h22771d5_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path  "versions.yml"               , emit: versions

    script:
    """
    MethylDackel extract \\
        $options.args \\
        $fasta \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
