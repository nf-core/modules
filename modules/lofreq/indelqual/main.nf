include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LOFREQ_INDELQUAL {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4"
    } else {
        container "quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4"
    }

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    lofreq indelqual \\
        $options.args \\
        -f $fasta \\
        -o ${prefix}.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
