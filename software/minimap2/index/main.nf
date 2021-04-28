// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MINIMAP2_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['']) }

    conda     (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"

    input:
    path(fasta)

    output:
    path("*.mmi")        , emit: index
    path "*.version.txt" ,emit: version

    script:
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta}.mmi \\
        $fasta
    ps
    minimap2 --version &> minimap2.version.txt
    """
}
