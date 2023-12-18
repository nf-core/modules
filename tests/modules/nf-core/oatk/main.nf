#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OATK } from '../../../../modules/nf-core/oatk/main.nf'

workflow test_oatk_mito {

    data_reads = Channel.of([[id:"ilDeiPorc1"],file(params.test_data['deilephila_porcellus']['mito']['hifi_reads'], checkIfExists: true)])
// TODO: update all paths to test_datasets    

//    hmm = file(params.test_data['deilephila_porcellus']['insecta_mito_hmm']['hifi_reads', checkIfExists: true)
    hmm = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/insecta_mito.fam")
    hmm_h3f = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/insecta_mito.fam.h3f")
    hmm_h3i = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/insecta_mito.fam.h3i")
    hmm_h3m = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/insecta_mito.fam.h3m")
    hmm_h3p = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/insecta_mito.fam.h3p")

// TODO: add nhmmscan to oatk conda/docker container
    nhmmscan = file("/lustre/scratch124/tol/projects/darwin/users/kk16/development/oatk/modules/tests/modules/nf-core/oatk/nhmmscan")
    OATK ( data_reads, hmm, hmm_h3f, hmm_h3i, hmm_h3m, hmm_h3p, [], [], [], [], [], nhmmscan, [])
}

workflow test_oatk_pltd {

    data_reads = Channel.of([[id:"ddAraThal4"],file("/lustre/scratch124/tol/projects/darwin/users/kk16/development/mitohifi/test/mito-datasets_original/plant-ddAraThal4-chloro/r/gbk.HiFiMapped.bam.filtered.first1000.fasta")])  

// TODO: update all paths to test_datasets
    hmm = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/embryophyta_pltd.fam")
    hmm_h3f = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/embryophyta_pltd.fam.h3f")
    hmm_h3i = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/embryophyta_pltd.fam.h3i")
    hmm_h3m = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/embryophyta_pltd.fam.h3m")
    hmm_h3p = file("/lustre/scratch124/tol/projects/darwin/users/cz3/organelle_asm/hmm_db/embryophyta_pltd.fam.h3p")

// TODO: add nhmmscan to oatk conda/docker container
    nhmmscan = file("/lustre/scratch124/tol/projects/darwin/users/kk16/development/oatk/modules/tests/modules/nf-core/oatk/nhmmscan")
    OATK ( data_reads, [], [], [], [], [], hmm, hmm_h3f, hmm_h3i, hmm_h3m, hmm_h3p, nhmmscan, [])
}
