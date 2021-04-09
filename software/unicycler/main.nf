// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNICYCLER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::unicycler=0.4.8' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3"
    } else {
        container "quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.scaffolds.fa'), emit: scaffolds
    tuple val(meta), path('*.assembly.gfa'), emit: gfa
    tuple val(meta), path('*.log')         , emit: log
    path  '*.version.txt'                  , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    unicycler \\
        --threads $task.cpus \\
        $options.args \\
        $input_reads \\
        --out ./

    mv assembly.fasta ${prefix}.scaffolds.fa
    mv assembly.gfa ${prefix}.assembly.gfa
    mv unicycler.log ${prefix}.unicycler.log

    echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//' > ${software}.version.txt
    """
}
