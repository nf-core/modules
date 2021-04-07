#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_COUNT } from '../../../../software/kallistobustools/count/main.nf' addParams( options: [args: "-x 10XV2 --workflow standard --cellranger --h5ad"] )
include { KALLISTOBUSTOOLS_COUNT as KALLISTOBUSTOOLS_COUNT_KITE} from '../../../../software/kallistobustools/count/main.nf' addParams( options: [args: "-x 10XV3 --workflow kite --cellranger --h5ad"] )

workflow test_kallistobustools_count {
    
    input   = [ [id:'test_standard'], // meta map 
                [file("${launchDir}/tests/data/delete_me/kallistobustools/count/SRR8599150_S1_L001_R1_001.20k_lines.fastq.gz", checkIfExists: true),
                 file("${launchDir}/tests/data/delete_me/kallistobustools/count/SRR8599150_S1_L001_R2_001.20k_lines.fastq.gz", checkIfExists: true)] 
                 ]  
    
    index   =   file("${launchDir}/tests/data/delete_me/kallistobustools/count/kb_ref.idx", checkIfExists: true)
    t2g     =   file("${launchDir}/tests/data/delete_me/kallistobustools/count/t2g.txt", checkIfExists: true)
    t1c     =   file('t1c_dummy')
    t2c     =   file('t2c_dummy')
    use_t1c =   false
    use_t2c =   false

    KALLISTOBUSTOOLS_COUNT (input,index,t2g,t1c,t2c,use_t1c,use_t2c)
}

workflow test_kallistobustools_count_kite {
    
    input   = [ [id:'test_lamanno'], // meta map
                [file("${launchDir}/tests/data/delete_me/kallistobustools/count/pbmc_1k_protein_v3_antibody_S2_L001_R1_001.100k_sub.fastq.gz", checkIfExists: true),
                 file("${launchDir}/tests/data/delete_me/kallistobustools/count/pbmc_1k_protein_v3_antibody_S2_L001_R2_001.100k_sub.fastq.gz", checkIfExists: true)] 
                 ]  
    
    index   =   file("${launchDir}/tests/data/delete_me/kallistobustools/count/mismatch.idx", checkIfExists: true)
    t2g     =   file("${launchDir}/tests/data/delete_me/kallistobustools/count/t2g.kite.txt", checkIfExists: true)
    t1c     =   file('t1c_dummy')
    t2c     =   file('t2c_dummy')
    use_t1c =   false
    use_t2c =   false

    KALLISTOBUSTOOLS_COUNT_KITE (input,index,t2g,t1c,t2c,use_t1c,use_t2c)
}
