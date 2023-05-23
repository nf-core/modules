process HLALA_TYPING {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hla-la=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'biocontainers/hla-la:1.0.3--hd03093a_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(graph)

    output:
    tuple val(meta), path("${meta.id}")                             , emit: results
    tuple val(meta), path("${meta.id}/extraction.bam*")             , emit: extraction
    tuple val(meta), path("${meta.id}/extraction_mapped.bam")       , emit: extraction_mapped
    tuple val(meta), path("${meta.id}/extraction_unmapped.bam")     , emit: extraction_unmpapped
    tuple val(meta), path("${meta.id}/hla/*")                       , emit: hla
    tuple val(meta), path("${meta.id}/*.fastq")                     , emit: fastq
    tuple val(meta), path("${meta.id}/reads_per_level.txt")         , emit: reads_per_level
    tuple val(meta), path("${meta.id}/remapped_with_a.bam*")        , emit: remapped

    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir graphs
    mv $graph graphs
    mv graphs /usr/local/opt/hla-la/

    mkdir $prefix

    /usr/local/opt/hla-la/src/HLA-LA.pl \\
        --BAM $bam \\
        --graph ../graphs/$graph \\
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
    mkdir ${prefix}/hla
    echo stub > ${prefix}/extraction.bam
    echo stub > ${prefix}/extraction.bam.bai
    echo stub > ${prefix}/extraction_mapped.bam
    echo stub > ${prefix}/extraction_unmapped.bam
    echo stub > ${prefix}/remapped_with_a.bam
    echo stub > ${prefix}/remapped_with_a.bam.bai
    echo stub > ${prefix}/R_1.fastq
    echo stub > ${prefix}/R_2.fastq
    echo stub > ${prefix}/R_U.fastq
    echo stub > ${prefix}/reads_per_level.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DPA1.txt
    echo stub > ${prefix}/hla/R1_pileup_G.txt
    echo stub > ${prefix}/hla/R1_pileup_E.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DRB3.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DRB1.txt
    echo stub > ${prefix}/hla/R1_readIDs_F.txt
    echo stub > ${prefix}/hla/R1_PP_DPA1_pairs.txt
    echo stub > ${prefix}/hla/R1_PP_DRB4_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_C.txt
    echo stub > ${prefix}/hla/R1_parameters.txt
    echo stub > ${prefix}/hla/R1_readIDs_G.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_V.txt
    echo stub > ${prefix}/hla/R1_readIDs_DRA.txt
    echo stub > ${prefix}/hla/R1_pileup_DRB4.txt
    echo stub > ${prefix}/hla/R1_pileup_DQB1.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_C.txt
    echo stub > ${prefix}/hla/R1_pileup_DRB3.txt
    echo stub > ${prefix}/hla/R1_readIDs_DRB4.txt
    echo stub > ${prefix}/hla/R1_PP_F_pairs.txt
    echo stub > ${prefix}/hla/R1_PP_E_pairs.txt
    echo stub > ${prefix}/hla/R1_readIDs_DQB1.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_A.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_H.txt
    echo stub > ${prefix}/hla/R1_readIDs_V.txt
    echo stub > ${prefix}/hla/R1_PP_DRB3_pairs.txt
    echo stub > ${prefix}/hla/R1_PP_A_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_DPA1.txt
    echo stub > ${prefix}/hla/R1_PP_G_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_K.txt
    echo stub > ${prefix}/hla/R1_PP_K_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_F.txt
    echo stub > ${prefix}/hla/R1_readIDs_K.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_K.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DRB4.txt
    echo stub > ${prefix}/hla/R1_pileup_B.txt
    echo stub > ${prefix}/hla/R1_pileup_V.txt
    echo stub > ${prefix}/hla/R1_pileup_H.txt
    echo stub > ${prefix}/hla/R1_readIDs_E.txt
    echo stub > ${prefix}/hla/R1_PP_DPB1_pairs.txt
    echo stub > ${prefix}/hla/summaryStatistics.txt
    echo stub > ${prefix}/hla/R1_readIDs_C.txt
    echo stub > ${prefix}/hla/R1_PP_H_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_DPB1.txt
    echo stub > ${prefix}/hla/R1_readIDs_DQA1.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DPB1.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_B.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DQB1.txt
    echo stub > ${prefix}/hla/R1_PP_V_pairs.txt
    echo stub > ${prefix}/hla/R1_readIDs_H.txt
    echo stub > ${prefix}/hla/R1_PP_DRB1_pairs.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_E.txt
    echo stub > ${prefix}/hla/R1_bestguess.txt
    echo stub > ${prefix}/hla/R1_PP_DRA_pairs.txt
    echo stub > ${prefix}/hla/R1_readIDs_DRB1.txt
    echo stub > ${prefix}/hla/R1_bestguess_G.txt
    echo stub > ${prefix}/hla/R1_PP_DQB1_pairs.txt
    echo stub > ${prefix}/hla/R1_readIDs_DRB3.txt
    echo stub > ${prefix}/hla/R1_readIDs_A.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DRA.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_G.txt
    echo stub > ${prefix}/hla/R1_pileup_DRA.txt
    echo stub > ${prefix}/hla/R1_readIDs_DPA1.txt
    echo stub > ${prefix}/hla/R1_PP_B_pairs.txt
    echo stub > ${prefix}/hla/R1_PP_DQA1_pairs.txt
    echo stub > ${prefix}/hla/R1_pileup_A.txt
    echo stub > ${prefix}/hla/R1_PP_C_pairs.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_F.txt
    echo stub > ${prefix}/hla/R1_columnIncompatibilities_DQA1.txt
    echo stub > ${prefix}/hla/R1_pileup_DRB1.txt
    echo stub > ${prefix}/hla/R1_readIDs_B.txt
    echo stub > ${prefix}/hla/histogram_matchesPerRead.txt
    echo stub > ${prefix}/hla/R1_pileup_DQA1.txt
    echo stub > ${prefix}/hla/R1_readIDs_DPB1.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hla-la: 1.0.3
    END_VERSIONS
    """
}
