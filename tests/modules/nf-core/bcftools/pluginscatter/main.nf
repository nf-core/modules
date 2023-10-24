#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_PLUGINSCATTER } from '../../../../../modules/nf-core/bcftools/pluginscatter/main.nf'

workflow test_bcftools_pluginscatter_sites_per_chunk {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ]

    BCFTOOLS_PLUGINSCATTER (
        input,
        100,
        [],
        [],
        [],
        []
    )
}

workflow test_bcftools_pluginscatter_scatter {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ]

    regions = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    BCFTOOLS_PLUGINSCATTER (
        input,
        [],
        "chr21",
        [],
        regions,
        []
    )
}

workflow test_bcftools_pluginscatter_scatter_file {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ]

    targets = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    scatter_file = Channel.of("chr21:6000000-41743940\tfile1", "chr21:41743941-46661900\tfile2")
        .collectFile(name:"scatter.tsv", newLine:true)

    BCFTOOLS_PLUGINSCATTER (
        input,
        [],
        [],
        scatter_file,
        [],
        targets
    )
}
