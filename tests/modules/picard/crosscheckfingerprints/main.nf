#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PICARD_CROSSCHECKFINGERPRINTS } from '../../../../modules/picard/crosscheckfingerprints/main.nf'

process STUB_HAPLOTYPE_MAP {

    output:
    path "haplotype_map.txt" ,emit:haplotype_map

    stub:
    """
    echo haplotype_map > haplotype_map.txt
    """
}

workflow test_picard_crosscheckfingerprints {

    input = [
        [ id:'test', single_end:false ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true), file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)],
    ]
    STUB_HAPLOTYPE_MAP ()
    PICARD_CROSSCHECKFINGERPRINTS ( input,[],STUB_HAPLOTYPE_MAP.out.haplotype_map )
}
