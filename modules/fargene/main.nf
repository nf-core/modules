// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.1'

process FARGENE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::fargene=0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fargene:0.1--py27h21c881e_4"
    } else {
        container "quay.io/biocontainers/fargene:0.1--py27h21c881e_4"
    }

    input:
    // input may be fasta (for genomes or longer contigs) or paired-end fastq (for metagenome), the latter in addition with --meta flag
    tuple val(meta), path(input)
    val hmm_model

    output:
    path "*.log"                                                                                 , emit: log
    path "${prefix}/results_summary.txt"                                                         , emit: txt
    tuple val(meta), path("${prefix}/hmmsearchresults/*.out")                                    , optional: true, emit: hmm
    tuple val(meta), path("${prefix}/predictedGenes/predicted-orfs.fasta")                       , optional: true, emit: orfs
    tuple val(meta), path("${prefix}/predictedGenes/predicted-orfs-amino.fasta")                 , optional: true, emit: orfs_amino
    tuple val(meta), path("${prefix}/predictedGenes/retrieved-contigs.fasta")                    , optional: true, emit: contigs
    tuple val(meta), path("${prefix}/predictedGenes/retrieved-contigs-peptides.fasta")           , optional: true, emit: contigs_pept
    tuple val(meta), path("${prefix}/predictedGenes/*filtered.fasta")                            , optional: true, emit: filtered
    tuple val(meta), path("${prefix}/predictedGenes/*filtered-peptides.fasta")                   , optional: true, emit: filtered_pept
    tuple val(meta), path("${prefix}/retrievedFragments/all_retrieved_*.fastq")                  , optional: true, emit: fragments
    tuple val(meta), path("${prefix}/retrievedFragments/retrievedFragments/trimmedReads/*.fasta"), optional: true, emit: trimmed
    tuple val(meta), path("${prefix}/spades_assembly/*")                                         , optional: true, emit: spades
    tuple val(meta), path("${prefix}/tmpdir/*.fasta")                                            , optional: true, emit: metagenome
    tuple val(meta), path("${prefix}/tmpdir/*.out")                                              , optional: true, emit: tmp
    path "versions.yml"                                                                          , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gzip \\
        -cdf $input \\
        > unziped.fa |
        fargene \\
            $options.args \\
            -p $task.cpus \\
            -i unziped.fa \\
            --hmm-model $hmm_model \\
            -o $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
