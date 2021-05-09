// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_SAMPLE {
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
    val sample_size

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "*.version.txt"               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        seqtk \\
            sample \\
            $options.args \\
            $reads \\
            $sample_size \\
            | gzip > ${prefix}.fastq.gz \\

        echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//' > ${software}.version.txt
        """
    } else {
        """
        seqtk \\
            sample \\
            $options.args \\
            $reads[0] \\
            $sample_size \\
            | gzip > ${prefix}_1.fastq.gz \\

        seqtk \\
            sample \\
            $options.args \\
            $reads[1] \\
            $sample_size \\
            | gzip > ${prefix}_2.fastq.gz \\

        echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//' > ${software}.version.txt
        """
    }
}
