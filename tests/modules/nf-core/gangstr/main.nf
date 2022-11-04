#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GANGSTR } from '../../../../modules/nf-core/gangstr/main.nf'

workflow test_gangstr_cram {

    bed = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
                .map({ bed ->
                    def content = bed.text.replace("\n","").tokenize("\t")
                    content.addAll(["5","CGCGC"])
                    [ content.join("\t") ]
                })
                .collectFile( newLine: true ) { content ->
                    [ "genome.bed", "${content[0]}" ]
                }

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]).combine(bed)

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    GANGSTR (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_gangstr_single_bam {

    bed = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
                .map({ bed ->
                    def content = bed.text.replace("\n","").tokenize("\t")
                    content.addAll(["5","CGCGC"])
                    [ content.join("\t") ]
                })
                .collectFile( newLine: true ) { content ->
                    [ "genome.bed", "${content[0]}" ]
                }

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]).combine(bed)

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    GANGSTR (
        input,
        fasta,
        fasta_fai
    )
}

workflow test_gangstr_multiple_bams {

    bed = Channel.value(file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true))
                .map({ bed ->
                    def content = bed.text.replace("\n","").tokenize("\t")
                    content.addAll(["5","CGCGC"])
                    [ content.join("\t") ]
                })
                .collectFile( newLine: true ) { content ->
                    [ "genome.bed", "${content[0]}" ]
                }

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true)
        ],
        [
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
        ]
    ]).combine(bed)

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    GANGSTR (
        input,
        fasta,
        fasta_fai
    )
}
