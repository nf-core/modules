// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QCAT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda     (params.enable_conda ? "bioconda::qcat=1.1.0" : null)
    container "quay.io/biocontainers/qcat:1.1.0--py_0"

    input:
    tuple val(meta), path(input_path)
    val(barcode_kit)
    
    output:
    tuple val(meta), path("fastq/*.fastq.gz") , emit: fastq
    path "*.version.txt"                      , emit: version

    script:
    """
    ## Unzip fastq file
    ## qcat doesnt support zipped files yet
    FILE=$input_path
    if [[ \$FILE == *.gz ]]
    then
    zcat $input_path > unzipped.fastq
    FILE=unzipped.fastq
    fi
    qcat  \\
    -f \$FILE \\
    -b ./fastq \\
    --kit $barcode_kit
            
    ## Zip fastq files (cannot find pigz command)
    gzip fastq/*
    qcat --version &> qcat.version.txt
    """
}
