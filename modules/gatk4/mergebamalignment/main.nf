// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MERGEBAMALIGNMENT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(aligned)
    path  unmapped
    path  fasta
    path  dict

    output:
    tuple val(meta), path('*.bam'), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gatk MergeBamAlignment \\
        ALIGNED=$aligned \\
        UNMAPPED=$unmapped \\
        R=$fasta \\
        O=${prefix}.bam \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
