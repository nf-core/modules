// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_MERGEPE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
    } else {
        container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        ln -s ${reads} ${prefix}.fastq.gz
        
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        seqtk \\
            mergepe \\
            $options.args \\
            ${reads} \\
            | gzip -n >> ${prefix}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
