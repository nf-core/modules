#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {VGAN_HAPLOCART} from '../../../../../modules/nf-core/vgan/haplocart/main.nf'

process STUB_HCFILES {

 output:
     path "hcfiles", emit: hcfiles

 script:
    """
    mkdir hcfiles
    touch hcfiles/graph.gg
    touch hcfiles/graph.og
    touch hcfiles/graph.xg
    touch hcfiles/graph.giraffe.gbz
    touch hcfiles/graph.dist
    touch hcfiles/graph_paths
    touch hcfiles/path_supports
    touch hcfiles/graph.gbwt
    touch hcfiles/children.txt
    touch hcfiles/parents.txt
    touch hcfiles/parsed_pangenome_mapping
    touch hcfiles/mappability.tsv
    touch hcfiles/k17_w18.min
    touch hcfiles/k31_w11.min
    touch versions.yml
    """
}

workflow test_vgan_haplocart_interleaved {

    STUB_HCFILES()

    input1 = [
        [ id:'test', single_end:true, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true)
             ]

    input2 = [ id:'test', single_end:true, format:'fastq'] // meta map

    VGAN_HAPLOCART(input1, input2, STUB_HCFILES.out.hcfiles)
}

workflow test_vgan_haplocart_paired_end_separate {

    STUB_HCFILES()
    rCRS_reads = file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true)
    rCRS_reads.copyTo(workDir + 'reads2.fq.gz')
    input1 = [
        [ id:'test', single_end:false, format:'fastq', reads2:file(workDir + "reads2.fq.gz", checkIfExists: true)], // meta map
          rCRS_reads
            ]

    input2 = [ id:'test', single_end:false, format:'fastq', reads2:file(workDir + "reads2.fq.gz", checkIfExists: true)] // meta map

    VGAN_HAPLOCART(input1, input2, STUB_HCFILES.out.hcfiles)
}

workflow test_vgan_haplocart_single_end {

    STUB_HCFILES()
    input1 = [
        [ id:'test', single_end:true, format:'fastq'], // meta map
        file(params.test_data['homo_sapiens']['illumina']['rCRS_reads'], checkIfExists: true)
             ]

    input2 = [id:'test', single_end:true, format:'fastq'] // meta map

    VGAN_HAPLOCART(input1, input2, STUB_HCFILES.out.hcfiles)
}
