// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GRAPHMAP2_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:['']) }

    conda (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0"
    } else {
        container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"
    }

    input:
    path fasta

    output:
    path "*.gmidx"      , emit: index
    path "versions.yml" , emit: version

    script:
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -I \\
        $options.args \\
        -r $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
