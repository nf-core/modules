// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::tiddit=2.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tiddit:2.12.1--py38h1773678_0"
    } else {
        container "quay.io/biocontainers/tiddit:2.12.1--py38h1773678_0"
    }

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.vcf"),         emit: vcf
    tuple val(meta), path("*.ploidy.tab"),  emit: ploidy
    tuple val(meta), path("*.signals.tab"), emit: signals
    path  "*.version.txt",                  emit: version

    script:
    def software = getSoftwareName(task.process)
    def output = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    def reference = fasta == "dummy_file.txt" ? "--ref $fasta" : ""
    """
    tiddit \\
        --sv $options.args \\
        --bam $bam \\
        $reference \\
        -o $output

    echo \$(tiddit -h 2>&1) | sed 's/^.*Version: //; s/(.*\$//' > ${software}.version.txt
    """
}
