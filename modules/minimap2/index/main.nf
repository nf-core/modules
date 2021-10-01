// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:['']) }

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/minimap2:2.21--h5bf99c6_0"
    }

    input:
    path fasta

    output:
    path "*.mmi"        , emit: index
    path "versions.yml" , emit: version

    script:
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $options.args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
