#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTK_BAFTEST      } from '../../../../../modules/nf-core/svtk/baftest/main.nf'
include { SVTK_VCF2BED      } from '../../../../../modules/nf-core/svtk/vcf2bed/main.nf'
include { MANTA_GERMLINE    } from '../../../../../modules/nf-core/manta/germline/main.nf'

workflow test_svtk_baftest {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    MANTA_GERMLINE (
        input,
        fasta,
        fasta_fai
    )

    SVTK_VCF2BED (
        MANTA_GERMLINE.out.diploid_sv_vcf.combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by:0)
    )

    baf_files = Channel.of([
        [ id:'test' ],
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/baftest/test.baf.txt.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/baftest/test.baf.txt.gz.tbi", checkIfExists: true)
    ])

    batch = Channel.of("sample\tgroup\tbatch","test\ttest\ttest").collectFile(name: "batch.txt", newLine:true, sort:true)

    SVTK_BAFTEST (
        SVTK_VCF2BED.out.bed.combine(baf_files, by:0).combine(batch)
    )
}
