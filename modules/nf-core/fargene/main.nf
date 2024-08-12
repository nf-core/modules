process FARGENE {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fargene:0.1--py27h21c881e_4' :
        'biocontainers/fargene:0.1--py27h21c881e_4' }"

    input:
    // input may be fasta (for genomes or longer contigs) or paired-end fastq (for metagenome), the latter in addition with --meta flag
    tuple val(meta), path(input)
    val hmm_model

    output:
    tuple val(meta), path("*.log")                                                                               , emit: log
    tuple val(meta), path("${prefix}/results_summary.txt")                                                       , emit: txt
    tuple val(meta), path("${prefix}/hmmsearchresults/*.out")                                    , optional: true, emit: hmm
    tuple val(meta), path("${prefix}/hmmsearchresults/retrieved-*.out")                          , optional: true, emit: hmm_genes
    tuple val(meta), path("${prefix}/predictedGenes/predicted-orfs.fasta")                       , optional: true, emit: orfs
    tuple val(meta), path("${prefix}/predictedGenes/predicted-orfs-amino.fasta")                 , optional: true, emit: orfs_amino
    tuple val(meta), path("${prefix}/predictedGenes/retrieved-contigs.fasta")                    , optional: true, emit: contigs
    tuple val(meta), path("${prefix}/predictedGenes/retrieved-contigs-peptides.fasta")           , optional: true, emit: contigs_pept
    tuple val(meta), path("${prefix}/predictedGenes/*filtered.fasta")                            , optional: true, emit: filtered
    tuple val(meta), path("${prefix}/predictedGenes/*filtered-peptides.fasta")                   , optional: true, emit: filtered_pept
    tuple val(meta), path("${prefix}/retrievedFragments/all_retrieved_*.fastq")                  , optional: true, emit: fragments
    tuple val(meta), path("${prefix}/retrievedFragments/trimmedReads/*.fasta")                   , optional: true, emit: trimmed
    tuple val(meta), path("${prefix}/spades_assembly/*")                                         , optional: true, emit: spades
    tuple val(meta), path("${prefix}/tmpdir/*.fasta")                                            , optional: true, emit: metagenome
    tuple val(meta), path("${prefix}/tmpdir/*.out")                                              , optional: true, emit: tmp
    path "versions.yml"                                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.log
    mkdir -p ${prefix}/{hmmsearchresults,predictedGenes,retrievedFragments}
    mkdir -p ${prefix}/retrievedFragments/trimmedReads/

    touch ${prefix}/results_summary.txt
    touch ${prefix}/hmmsearchresults/retrieved-${prefix}.out
    touch ${prefix}/hmmsearchresults/${prefix}.out
    touch ${prefix}/predictedGenes/predicted-orfs.fasta
    touch ${prefix}/predictedGenes/predicted-orfs-amino.fasta
    touch ${prefix}/predictedGenes/retrieved-contigs.fasta
    touch ${prefix}/predictedGenes/retrieved-contigs-peptides.fasta
    touch ${prefix}/predictedGenes/${prefix}-filtered.fasta
    touch ${prefix}/predictedGenes/${prefix}-filtered-peptides.fasta
    touch ${prefix}/retrievedFragments/all_retrieved_${prefix}.fastq
    touch ${prefix}/retrievedFragments/trimmedReads/${prefix}.fasta


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fargene: $VERSION
    END_VERSIONS
    """

}
