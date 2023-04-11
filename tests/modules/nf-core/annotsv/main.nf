#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANNOTSV } from '../../../../modules/nf-core/annotsv/main.nf'

workflow test_annotsv {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data["homo_sapiens"]["illumina"]["test_sv_vcf"], checkIfExists: true),
        file(params.test_data["homo_sapiens"]["illumina"]["test_sv_vcf_tbi"], checkIfExists: true)
    ]

    // annotations = [
    //     [ id: 'annotations' ],
    //     [] // For stub use only, this will fail if the actual module is run like this
    // ]

    annotations = [
        [ id:'annotations' ],
        file("/home/nvnieuwk/Documents/data/AnnotSV/share/AnnotSV", checkIfExists: true)
    ]

    genes = Channel
        .of('GENE1', 'GENE2', 'GENE3')
        .collectFile(name:'gene_candidates.txt', newLine:true)
        .map { [[id:'test'], it]}

    small_variants = [
        [ id:'test' ],
        file(params.test_data["homo_sapiens"]["illumina"]["test2_haplotc_vcf_gz"], checkIfExists:true)
    ]

    false_positives = [
        [ id:'test' ],
        file(params.test_data["homo_sapiens"]["illumina"]["test_haplotc_cnn_vcf_gz"], checkIfExists:true)
    ]

    gene_transcripts = Channel
        .of('GENE1 GENE2 GENE3')
        .collectFile(name:'gene_transcripts.txt')
        .map { [[id:'test'], it]}

    ANNOTSV (
        input,
        annotations,
        genes,
        small_variants,
        false_positives,
        gene_transcripts
    )
}
