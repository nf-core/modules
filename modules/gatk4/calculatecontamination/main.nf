process GATK4_CALCULATECONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(pileup), path(matched)
    val segmentout

    output:
    tuple val(meta), path('*.contamination.table')               , emit: contamination
    tuple val(meta), path('*.segmentation.table') , optional:true, emit: segmentation
    path "versions.yml"                                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def matched_command = matched ? " -matched ${matched} " : ''
    def segment_command = segmentout ? " -segments ${prefix}.segmentation.table" : ''
    """
    gatk CalculateContamination \\
        -I $pileup \\
        $matched_command \\
        -O ${prefix}.contamination.table \\
        $segment_command \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
