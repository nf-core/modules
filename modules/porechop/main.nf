// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PORECHOP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::porechop=0.2.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2"
    } else {
        container "quay.io/biocontainers/porechop:0.2.4--py38h8c62d01_2"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_porechop.fastq.gz")  , emit: reads
    path "versions.yml"                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    porechop \\
        -i ${reads} \\
        -t ${task.cpus} \\
        ${options.args} \\
        -o ${prefix}_porechop.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( porechop --version )
    END_VERSIONS
    """
}
