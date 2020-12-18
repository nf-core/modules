// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hed8969a_1'
    } else {
        container 'quay.io/biocontainers/bowtie:1.3.0--py38hed8969a_1'
    }

    input:
    path fasta

    output:
    path 'bowtie', emit: index
    path '*.version.txt', emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    mkdir bowtie
    bowtie-build --threads $task.cpus $fasta bowtie/${fasta.baseName}
    echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
