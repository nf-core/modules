def VERSION = '0.1' // Version information not provided by tool on CLI

process FARGENE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fargene=0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fargene:0.1--py27h21c881e_4' :
        'quay.io/biocontainers/fargene:0.1--py27h21c881e_4' }"

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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    fargene \\
        $args \\
        -p $task.cpus \\
        -i $input \\
        --hmm-model $hmm_model \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fargene: $VERSION
    END_VERSIONS
    """
}
