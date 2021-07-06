// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MERGEVCFS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::gatk4=4.1.9.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0'
    } else {
        container 'quay.io/biocontainers/gatk4:4.1.9.0--py39_0'
    }

    input:
    tuple val(meta), path(vcfs)
    path  ref_dict
    val   use_ref_dict

    output:
    tuple val(meta), path('*.vcf.gz'), emit: vcf
    path  '*.version.txt'            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Make list of VCFs to merge
    def input = ""
    for (vcf in vcfs) {
        input += " I=${vcf}"
    }
    def ref = use_ref_dict ? "D=${ref_dict}" : ""
    """
    gatk MergeVcfs \\
        $input \\
        O=${prefix}.vcf.gz \\
        $ref \\
        $options.args

    echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > ${software}.version.txt
    """
}
