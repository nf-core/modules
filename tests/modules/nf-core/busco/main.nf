#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../../modules/nf-core/busco/main.nf'

workflow test_busco_genome_single_fasta {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10', // Launch with 'auto' to use --auto-lineage, and specified lineages // 'auto' removed from test due to memory issues
        [], // Download busco lineage
        [] // No config
    )

}

workflow test_busco_genome_multi_fasta {

    input = [
        [ id:'test' ], // meta map
        [
            file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
            file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

}

workflow test_busco_eukaryote_metaeuk {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [] // No config
    )

}

workflow test_busco_eukaryote_augustus {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [] // No config
    )

}

workflow test_busco_protein {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

}

workflow test_busco_transcriptome {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

}
