// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PRODIGAL {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::prodigal=2.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/prodigal:2.6.3--h516909a_2"
    } else {
        container "quay.io/biocontainers/prodigal:2.6.3--h516909a_2"
    }

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}"), emit: gene_annotations
    tuple val(meta), path("${prefix}.fna"), emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa"), emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt"), emit: all_gene_annotations
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    prodigal -i "${genome}" \\
        $options.args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """
}
