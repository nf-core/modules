#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VGAN_HAPLOCART } from '../../../../../modules/nf-core/vgan/haplocart/main.nf' addParams( options: [args: ''] )

    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.gg", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.og", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.xg", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.giraffe.gbz", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.dist", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph_paths", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/path_supports", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.gbwt", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/children.txt", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/parents.txt", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/parsed_pangenome_mapping", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/mappability.tsv", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/k17_w18.min", checkIfExists: true)
    file("ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/k31_w11.min", checkIfExists: true)


workflow test_vgan_haplocart_paired_end_interleaved {

    input = [
        [ id:'test', single_end:false, format:'fastq'], // meta map
        file('hcfiles', checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        []
            ]

    VGAN_HAPLOCART(input)
}

workflow test_vgan_haplocart_paired_end_separate {

    rCRS_reads = file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true)
    rCRS_reads.copyTo(workDir + 'reads2.fq.gz')
    input = [
        [ id:'test', single_end:false, format:'fastq'], // meta map
        file('hcfiles', checkIfExists: true),
        rCRS_reads,
        file(workDir + "reads2.fq.gz", checkIfExists: true)
            ]

    VGAN_HAPLOCART(input)
}

workflow test_vgan_haplocart_single_end {

    input = [
        [ id:'test', single_end:true, format:'fastq'], // meta map
        file('hcfiles', checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true),
        []
        ]

    VGAN_HAPLOCART(input)
}
