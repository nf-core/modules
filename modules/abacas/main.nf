// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ABACAS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::abacas=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abacas:1.3.1--pl526_0"
    } else {
        container "quay.io/biocontainers/abacas:1.3.1--pl526_0"
    }

    input:
    tuple val(meta), path(scaffold)
    path  fasta

    output:
    tuple val(meta), path('*.abacas*'), emit: results
    path '*.version.txt'              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    abacas.pl \\
        -r $fasta \\
        -q $scaffold \\
        $options.args \\
        -o ${prefix}.abacas

    mv nucmer.delta ${prefix}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${prefix}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${prefix}.abacas.nucmer.tiling
    mv unused_contigs.out ${prefix}.abacas.unused.contigs.out
    echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//' > ${software}.version.txt
    """
}
