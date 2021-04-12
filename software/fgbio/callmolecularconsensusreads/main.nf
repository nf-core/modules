include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FGBIO_CALLMOLECULARCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fgbio:1.3.0--0"
    } else {
        container "quay.io/biocontainers/fgbio:1.3.0--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    fgbio \\
        CallMolecularConsensusReads \\
        -i $bam \\
        $options.args \\
        -o ${prefix}.bam
    fgbio --version | sed -e "s/fgbio v//g" > ${software}.version.txt
    """
}
