// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::bowtie=1.3.0" : null)
    container "quay.io/biocontainers/bowtie:1.3.0--py38hed8969a_1"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    
    output:
    tuple val(meta), path("*.sam")            , emit: sam
    tuple val(meta), path("*.out")            , emit: log
    path  "bowtie.version.txt"                , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    index_array = index.collect()
    index = index_array[0].baseName - ~/.\d$/
    """
    bowtie $options.args  \\
    --threads $task.cpus \\
    ${index}  \\
    -q ${reads} \\
    --un ${prefix}.unAl > ${prefix}.sam 2> ${prefix}.out
    bowtie --version | head -n 1 | cut -d" " -f3 > ${software}.version.txt
    """
}