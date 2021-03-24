// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSOR_MSI {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2"
    } else {
        container "quay.io/biocontainers/msisensor:0.5--hb3646a4_2"
    }

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai), val(metascan), path(homopolymers)

    output:
    tuple val(meta), path("*.{output,output_dis,output_germline,output_somatic}"), emit: report

    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    msisensor msi -d ${homopolymers} \\
                    -n ${normal_bam} \\
                    -t ${tumor_bam} \\
                    -o ${prefix}.paired.output \\
                    $options.args

    msisensor 2>&1 |sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p' > ${software}.version.txt
    """
}
