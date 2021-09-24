// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process QCAT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::qcat=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/qcat:1.1.0--py_0"
    } else {
        container "quay.io/biocontainers/qcat:1.1.0--py_0"
    }

    input:
    tuple val(meta), path(reads)
    val   barcode_kit

    output:
    tuple val(meta), path("fastq/*.fastq.gz"), emit: reads
    path "versions.yml"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ## Unzip fastq file
    ## qcat doesn't support zipped files yet
    FILE=$reads
    if [[ \$FILE == *.gz ]]
    then
        zcat $reads > unzipped.fastq
        FILE=unzipped.fastq
    fi

    qcat \\
        -f \$FILE \\
        -b ./fastq \\
        --kit $barcode_kit

    ## Zip fastq files
    gzip fastq/*

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(qcat --version 2>&1 | sed 's/^.*qcat //; s/ .*\$//')
    END_VERSIONS
    """
}
