// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_CALCULATECONTAMINATION {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(pileup), path(matched)

    val segmentout

    output:
    tuple val(meta), path('*.contamination.table'), emit: contamination
    tuple val(meta), path('*.segmentation.table'), optional:true, emit: segmentation
    path "versions.yml"           , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def matchedCommand = ''
    def segmentCommand = ''

    matchedCommand = matched ? " -matched ${matched} " : ''

    if(segmentout){
        segmentCommand = " -segments ${prefix}.segmentation.table"
    } else {
        segmentCommand = ''
    }

    """
    gatk CalculateContamination \\
        -I $pileup \\
        ${matchedCommand} \\
        -O ${prefix}.contamination.table \\
        ${segmentCommand} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(gatk --version 2>&1 | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
