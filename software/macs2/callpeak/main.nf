// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process MACS2_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/macs2:2.2.7.1--py37h516909a_0"
    //container "https://depot.galaxyproject.org/singularity/macs2:2.2.7.1--py37h516909a_0"

    conda (params.conda ? "bioconda::macs2=2.2.7.1" : null)

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val macs2_gsize
    val options

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls"), emit: xls
    tuple val(meta), path("*.gappedPeak"), emit: gapped optional true
    tuple val(meta), path("*.bed"), emit: bed optional true
    tuple val(meta), path("*.bdg"), emit: bdg optional true
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def format   = meta.single_end ? 'BAM' : 'BAMPE'
    def control  = controlbam ? "--control $controlbam" : ''
    """
    macs2 \\
        callpeak \\
        $ioptions.args \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $ipbam \\
         $control

    macs2 --version | sed -e "s/macs2 //g" > ${software}.version.txt
    """
}
