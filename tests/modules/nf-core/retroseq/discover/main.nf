#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RETROSEQ_DISCOVER } from '../../../../../modules/nf-core/retroseq/discover/main.nf'

// Add a process to generate elements.tab
process generate_elements_tab {
    output:
    path 'elements.txt'
    script:
    """
    echo -e 'ALU'\t'annotated.bed' > elements.txt
    """
}

workflow test_retroseq_discover {
    fileee = generate_elements_tab()

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    tab = [
        [ id:'elementlist'], //meta map
        file(params.test_data['homo_sapiens']['genome']['genome_21_annotated_bed'], checkIfExists: true)
    ]

    tablist = [
            [ id:'tablist'], //meta map
            file(generate_elements_tab.view)
    ]

    // Connect the process to the workflow
    RETROSEQ_DISCOVER(bam, tablist)
}

workflow test_retroseq_discover_stub {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    ]

    RETROSEQ_DISCOVER (bam, tablist)
}
