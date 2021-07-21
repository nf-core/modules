#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KALLISTOBUSTOOLS_COUNT } from '../../../../modules/kallistobustools/count/main.nf' addParams( options: [args:"--cellranger"] )

workflow test_kallistobustools_count {
    
    input   = [ [id:'test_standard'], // meta map 
                [file("https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/10xgenomics/test_1.fastq.gz?raw=true", checkIfExists: true),
                 file("https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/10xgenomics/test_2.fastq.gz?raw=true", checkIfExists: true)] 
                 ]  
    
    index       =   file("https://github.com/FloWuenne/test-datasets/blob/scrnaseq/reference/kallistobustools/kb_ref.idx?raw=true", checkIfExists: true)
    t2g         =   file("https://raw.githubusercontent.com/FloWuenne/test-datasets/scrnaseq/reference/kallistobustools/t2g.txt", checkIfExists: true)
    t1c         =   file('t1c_dummy')
    t2c         =   file('t2c_dummy')
    use_t1c     =   false
    use_t2c     =   false
    workflow    =   "standard"
    technology  =   "10XV3"

    KALLISTOBUSTOOLS_COUNT (input,index,t2g,t1c,t2c,use_t1c,use_t2c,workflow,technology)
}
