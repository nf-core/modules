#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf'

workflow test_busco_genome_lineage {
    input = [ [ id:'test', mode:"genome", lineage:"bacteroidales_odb10", autolineage:""],
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    BUSCO(input)
}

workflow test_busco_genome_autolineage {
    input = [ [ id:'test', mode:"genome", lineage:"", autolineage:"auto"],
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    BUSCO(input)
}

workflow test_busco_genome {
    input = [ [ id:'test', mode:"genome", lineage:"", autolineage:""],
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    BUSCO(input)
}

/*workflow test_busco_protein {
    input = [ [ id:'test', mode:"protein", lineage:"", autolineage:""],
        file(params.test_data['bacteroides_fragilis']['illumina']['test1.contigs.fa.gz'], checkIfExists: true) ]

    BUSCO(input)
}*/

workflow test_busco_transcriptome {
    input = [ [ id:'test', mode:"transcriptome", lineage:"", autolineage:""],
        file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true) ]

    BUSCO(input)
}
