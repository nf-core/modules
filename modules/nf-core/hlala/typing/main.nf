process HLALA_TYPING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.4--h077b44d_1':
        'biocontainers/hla-la:1.0.4--h077b44d_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(graph)

    output:
    tuple val(meta), path("results")                             , emit: results
    tuple val(meta), path("results/extraction.bam")              , emit: extraction
    tuple val(meta), path("results/extraction.bam.bai")          , emit: extraction_index
    tuple val(meta), path("results/extraction_mapped.bam")       , emit: extraction_mapped
    tuple val(meta), path("results/extraction_unmapped.bam")     , emit: extraction_unmpapped
    tuple val(meta), path("results/hla/*")                       , emit: hla
    tuple val(meta), path("results/*.fastq")                     , emit: fastq
    tuple val(meta), path("results/reads_per_level.txt")         , emit: reads_per_level
    tuple val(meta), path("results/remapped_with_a.bam")         , emit: remapped
    tuple val(meta), path("results/remapped_with_a.bam.bai")     , emit: remapped_index
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def bin = ""
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        bin="\$CONDA_PREFIX/opt/hla-la/src/HLA-LA.pl"
    } else {
        bin="/usr/local/opt/hla-la/src/HLA-LA.pl"
    }

    """
    ${bin} \\
        --BAM $bam \\
        --customGraphDir ${graph} \\
        --sampleID $prefix \\
        --workingDir . \\
        --maxThreads $task.cpus \\
        $args

    mkdir -p results
    mv ${prefix}/ results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: ${VERSION}
    END_VERSIONS
    """

    stub:
    def VERSION = '1.0.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p results

    touch results/R_1.fastq
    touch results/R_2.fastq
    touch results/R_U.fastq
    touch results/extraction.bam
    touch results/extraction.bam.bai
    touch results/extraction_mapped.bam
    touch results/extraction_unmapped.bam
    touch results/remapped_with_a.bam
    touch results/remapped_with_a.bam.bai
    touch results/reads_per_level.txt

    mkdir results/hla
    touch results/hla/histogram_matchesPerRead.txt
    touch results/hla/R1_bestguess_G.txt
    touch results/hla/R1_bestguess.txt
    touch results/hla/R1_columnIncompatibilities_A.txt
    touch results/hla/R1_columnIncompatibilities_B.txt
    touch results/hla/R1_columnIncompatibilities_C.txt
    touch results/hla/R1_columnIncompatibilities_DPA1.txt
    touch results/hla/R1_columnIncompatibilities_DPB1.txt
    touch results/hla/R1_columnIncompatibilities_DQA1.txt
    touch results/hla/R1_columnIncompatibilities_DQB1.txt
    touch results/hla/R1_columnIncompatibilities_DRA.txt
    touch results/hla/R1_columnIncompatibilities_DRB1.txt
    touch results/hla/R1_columnIncompatibilities_DRB3.txt
    touch results/hla/R1_columnIncompatibilities_DRB4.txt
    touch results/hla/R1_columnIncompatibilities_E.txt
    touch results/hla/R1_columnIncompatibilities_F.txt
    touch results/hla/R1_columnIncompatibilities_G.txt
    touch results/hla/R1_columnIncompatibilities_H.txt
    touch results/hla/R1_columnIncompatibilities_K.txt
    touch results/hla/R1_columnIncompatibilities_V.txt
    touch results/hla/R1_parameters.txt
    touch results/hla/R1_pileup_A.txt
    touch results/hla/R1_pileup_B.txt
    touch results/hla/R1_pileup_C.txt
    touch results/hla/R1_pileup_DPA1.txt
    touch results/hla/R1_pileup_DPB1.txt
    touch results/hla/R1_pileup_DQA1.txt
    touch results/hla/R1_pileup_DQB1.txt
    touch results/hla/R1_pileup_DRA.txt
    touch results/hla/R1_pileup_DRB1.txt
    touch results/hla/R1_pileup_DRB3.txt
    touch results/hla/R1_pileup_DRB4.txt
    touch results/hla/R1_pileup_E.txt
    touch results/hla/R1_pileup_F.txt
    touch results/hla/R1_pileup_G.txt
    touch results/hla/R1_pileup_H.txt
    touch results/hla/R1_pileup_K.txt
    touch results/hla/R1_pileup_V.txt
    touch results/hla/R1_PP_A_pairs.txt
    touch results/hla/R1_PP_B_pairs.txt
    touch results/hla/R1_PP_C_pairs.txt
    touch results/hla/R1_PP_DPA1_pairs.txt
    touch results/hla/R1_PP_DPB1_pairs.txt
    touch results/hla/R1_PP_DQA1_pairs.txt
    touch results/hla/R1_PP_DQB1_pairs.txt
    touch results/hla/R1_PP_DRA_pairs.txt
    touch results/hla/R1_PP_DRB1_pairs.txt
    touch results/hla/R1_PP_DRB3_pairs.txt
    touch results/hla/R1_PP_DRB4_pairs.txt
    touch results/hla/R1_PP_E_pairs.txt
    touch results/hla/R1_PP_F_pairs.txt
    touch results/hla/R1_PP_G_pairs.txt
    touch results/hla/R1_PP_H_pairs.txt
    touch results/hla/R1_PP_K_pairs.txt
    touch results/hla/R1_PP_V_pairs.txt
    touch results/hla/R1_readIDs_A.txt
    touch results/hla/R1_readIDs_B.txt
    touch results/hla/R1_readIDs_C.txt
    touch results/hla/R1_readIDs_DPA1.txt
    touch results/hla/R1_readIDs_DPB1.txt
    touch results/hla/R1_readIDs_DQA1.txt
    touch results/hla/R1_readIDs_DQB1.txt
    touch results/hla/R1_readIDs_DRA.txt
    touch results/hla/R1_readIDs_DRB1.txt
    touch results/hla/R1_readIDs_DRB3.txt
    touch results/hla/R1_readIDs_DRB4.txt
    touch results/hla/R1_readIDs_E.txt
    touch results/hla/R1_readIDs_F.txt
    touch results/hla/R1_readIDs_G.txt
    touch results/hla/R1_readIDs_H.txt
    touch results/hla/R1_readIDs_K.txt
    touch results/hla/R1_readIDs_V.txt
    touch results/hla/summaryStatistics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: ${VERSION}
    END_VERSIONS
    """
}
