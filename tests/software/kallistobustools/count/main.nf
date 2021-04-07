#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_COUNT } from '../../../../software/kallistobustools/count/main.nf' addParams( options: [args: "-x 10XV2 --workflow standard --cellranger --h5ad"] )

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
