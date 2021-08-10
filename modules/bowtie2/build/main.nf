// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0'
    } else {
        container 'quay.io/biocontainers/bowtie2:2.4.4--py36hd4290be_0'
    }

    input:
    path fasta

    output:
    path 'bowtie2'      , emit: index
    path '*.version.txt', emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    mkdir bowtie2
    bowtie2-build $options.args --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
