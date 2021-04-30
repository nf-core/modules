// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GRAPHMAP2_INDEX {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"

    input:
    path(fasta)

    output:
    path("*.gmidx")       ,emit: index
    path "*.version.txt"  ,emit: version

    script:
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -I \\
        -r $fasta
    echo \$(graphmap2 2>&1) > graphmap2.version.txt
    """
}
