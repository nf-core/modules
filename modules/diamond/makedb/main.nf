// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DIAMOND_MAKEDB {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? 'bioconda::diamond=2.0.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container 'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0'
    }

    input:
    path fasta

    output:
    path "${fasta}.dmnd", emit: db
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    diamond makedb \\
        --threads $task.cpus \\
        --in  $fasta \\
        -d $fasta \\
        $options.args
    echo \$(diamond --version 2>&1) | tail -n 1 | sed 's/^diamond version //' > ${software}.version.txt
    """
}
