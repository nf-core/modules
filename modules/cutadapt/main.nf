// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1'
    } else {
        container 'quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    cutadapt \\
        --cores $task.cpus \\
        $options.args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cutadapt --version)
    END_VERSIONS
    """
}
