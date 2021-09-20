// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LIMA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::lima=2.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lima:2.2.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/lima:2.2.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(ccs)
    path primers

    output:
    tuple val(meta), path("*.bam")    , emit: bam
    tuple val(meta), path("*.bam.pbi"), emit: pbi
    tuple val(meta), path("*.xml")    , emit: xml
    tuple val(meta), path("*.json")   , emit: json
    tuple val(meta), path("*.clips")  , emit: clips
    tuple val(meta), path("*.counts") , emit: counts
    tuple val(meta), path("*.guess")  , emit: guess
    tuple val(meta), path("*.report") , emit: report
    tuple val(meta), path("*.summary"), emit: summary
    path("*.version.txt")             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def lima_out = ccs.toString().replaceAll(/bam$/, 'fl.bam')
    """
    lima \\
        $ccs \\
        $primers \\
        $lima_out \\
        -j $task.cpus \\
        $options.args

    lima --version | sed 's/lima //g' | sed 's/ (.\\+//g' > ${software}.version.txt
    """
}
