// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_ALN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8"
    } else {
        container "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path(reads), path("*.sai"), emit: sai
    path "*.version.txt"                       , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        bwa aln \\
            $options.args \\
            -t $task.cpus \\
            -f ${prefix}.sai \\
            \$INDEX \\
            ${reads}

        echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
        """
    } else {
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

        bwa aln \\
            $options.args \\
            -t $task.cpus \\
            -f ${prefix}.1.sai \\
            \$INDEX \\
            ${reads[0]}

        bwa aln \\
            $options.args \\
            -t $task.cpus \\
            -f ${prefix}.2.sai \\
            \$INDEX \\
            ${reads[1]}

        echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
        """
    }
}
