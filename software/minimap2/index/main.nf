// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:['']) }

    conda (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3"
    } else {
        container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    }

    input:
    path fasta

    output:
    path "*.mmi"        , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $options.args \\
        $fasta

    echo \$(minimap2 --version 2>&1) > ${software}.version.txt
    """
}
