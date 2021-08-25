include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ULTRA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ultra_bioinformatics=0.0.3.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.0.3.3--pyh5e36f6f_1"
    } else {
        container "quay.io/biocontainers/ultra_bioinformatics:0.0.3.3--pyh5e36f6f_1"
    }

    input:
    tuple val(meta), path(reads)
    path(genome)
    path(gtf)

    output:
    tuple val(meta), path("outfolder/*.sam"), emit: sam
    path "*.version.txt"                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    uLTRA pipeline \\
        --t $task.cpus \\
        --prefix $prefix \\
        $options.args \\
        \$(pwd)/$genome \\
        \$(pwd)/$gtf \\
        \$(pwd)/$reads \\
        outfolder/

    echo \$(uLTRA --version 2>&1) > ${software}.version.txt
    """
}
