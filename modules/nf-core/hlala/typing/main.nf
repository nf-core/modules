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
    tuple val(meta), path("${meta.id}")                             , emit: results
    tuple val(meta), path("${meta.id}/extraction.bam")              , emit: extraction
    tuple val(meta), path("${meta.id}/extraction.bam.bai")          , emit: extraction_index
    tuple val(meta), path("${meta.id}/extraction_mapped.bam")       , emit: extraction_mapped
    tuple val(meta), path("${meta.id}/extraction_unmapped.bam")     , emit: extraction_unmpapped
    tuple val(meta), path("${meta.id}/hla/*")                       , emit: hla
    tuple val(meta), path("${meta.id}/*.fastq")                     , emit: fastq
    tuple val(meta), path("${meta.id}/reads_per_level.txt")         , emit: reads_per_level
    tuple val(meta), path("${meta.id}/remapped_with_a.bam")         , emit: remapped
    tuple val(meta), path("${meta.id}/remapped_with_a.bam.bai")     , emit: remapped_index
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def bin = ""
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        bin="\$CONDA_PREFIX/opt/hla-la/src/HLA-LA.pl"
    } else {
        bin="/usr/local/opt/hla-la/src/HLA-LA.pl"
    }

    """
    ${bin} \\
        --BAM $bam \\
        --graph PRG_MHC_GRCh38_withIMGT \\
        --customGraphDir ${graph} \\
        --sampleID $prefix \\
        --workingDir . \\
        --maxThreads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: 1.0.3
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    touch ${prefix}/R_1.fastq
    touch ${prefix}/R_2.fastq
    touch ${prefix}/R_U.fastq
    touch ${prefix}/extraction.bam
    touch ${prefix}/extraction.bam.bai
    touch ${prefix}/extraction_mapped.bam
    touch ${prefix}/extraction_unmapped.bam
    touch ${prefix}/remapped_with_a.bam
    touch ${prefix}/remapped_with_a.bam.bai
    touch ${prefix}/reads_per_level.txt

    mkdir ${prefix}/hla
    touch ${prefix}/hla/histogram_matchesPerRead.txt
    touch ${prefix}/hla/R1_bestguess_G.txt
    touch ${prefix}/hla/R1_bestguess.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_A.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_B.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_C.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DPA1.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DPB1.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DQA1.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DQB1.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DRA.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DRB1.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DRB3.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_DRB4.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_E.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_F.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_G.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_H.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_K.txt
    touch ${prefix}/hla/R1_columnIncompatibilities_V.txt
    touch ${prefix}/hla/R1_parameters.txt
    touch ${prefix}/hla/R1_pileup_A.txt
    touch ${prefix}/hla/R1_pileup_B.txt
    touch ${prefix}/hla/R1_pileup_C.txt
    touch ${prefix}/hla/R1_pileup_DPA1.txt
    touch ${prefix}/hla/R1_pileup_DPB1.txt
    touch ${prefix}/hla/R1_pileup_DQA1.txt
    touch ${prefix}/hla/R1_pileup_DQB1.txt
    touch ${prefix}/hla/R1_pileup_DRA.txt
    touch ${prefix}/hla/R1_pileup_DRB1.txt
    touch ${prefix}/hla/R1_pileup_DRB3.txt
    touch ${prefix}/hla/R1_pileup_DRB4.txt
    touch ${prefix}/hla/R1_pileup_E.txt
    touch ${prefix}/hla/R1_pileup_F.txt
    touch ${prefix}/hla/R1_pileup_G.txt
    touch ${prefix}/hla/R1_pileup_H.txt
    touch ${prefix}/hla/R1_pileup_K.txt
    touch ${prefix}/hla/R1_pileup_V.txt
    touch ${prefix}/hla/R1_PP_A_pairs.txt
    touch ${prefix}/hla/R1_PP_B_pairs.txt
    touch ${prefix}/hla/R1_PP_C_pairs.txt
    touch ${prefix}/hla/R1_PP_DPA1_pairs.txt
    touch ${prefix}/hla/R1_PP_DPB1_pairs.txt
    touch ${prefix}/hla/R1_PP_DQA1_pairs.txt
    touch ${prefix}/hla/R1_PP_DQB1_pairs.txt
    touch ${prefix}/hla/R1_PP_DRA_pairs.txt
    touch ${prefix}/hla/R1_PP_DRB1_pairs.txt
    touch ${prefix}/hla/R1_PP_DRB3_pairs.txt
    touch ${prefix}/hla/R1_PP_DRB4_pairs.txt
    touch ${prefix}/hla/R1_PP_E_pairs.txt
    touch ${prefix}/hla/R1_PP_F_pairs.txt
    touch ${prefix}/hla/R1_PP_G_pairs.txt
    touch ${prefix}/hla/R1_PP_H_pairs.txt
    touch ${prefix}/hla/R1_PP_K_pairs.txt
    touch ${prefix}/hla/R1_PP_V_pairs.txt
    touch ${prefix}/hla/R1_readIDs_A.txt
    touch ${prefix}/hla/R1_readIDs_B.txt
    touch ${prefix}/hla/R1_readIDs_C.txt
    touch ${prefix}/hla/R1_readIDs_DPA1.txt
    touch ${prefix}/hla/R1_readIDs_DPB1.txt
    touch ${prefix}/hla/R1_readIDs_DQA1.txt
    touch ${prefix}/hla/R1_readIDs_DQB1.txt
    touch ${prefix}/hla/R1_readIDs_DRA.txt
    touch ${prefix}/hla/R1_readIDs_DRB1.txt
    touch ${prefix}/hla/R1_readIDs_DRB3.txt
    touch ${prefix}/hla/R1_readIDs_DRB4.txt
    touch ${prefix}/hla/R1_readIDs_E.txt
    touch ${prefix}/hla/R1_readIDs_F.txt
    touch ${prefix}/hla/R1_readIDs_G.txt
    touch ${prefix}/hla/R1_readIDs_H.txt
    touch ${prefix}/hla/R1_readIDs_K.txt
    touch ${prefix}/hla/R1_readIDs_V.txt
    touch ${prefix}/hla/summaryStatistics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: 1.0.3
    END_VERSIONS
    """
}
