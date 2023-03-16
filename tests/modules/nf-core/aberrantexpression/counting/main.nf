#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESSGENEANNOTATION  } from '../../../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { COUNTREADS                } from '../../../../../modules/nf-core/aberrantexpression/counting/countreads'
include { MERGECOUNTS               } from '../../../../../modules/nf-core/aberrantexpression/counting/mergecounts'
include { FILTERCOUNT               } from '../../../../../modules/nf-core/aberrantexpression/counting/filtercount'

// TODO-csandu: Add as separate file when creating the docker image
params.config='/usr/local/lib/python3.11/site-packages/drop/template/config.yaml'

process prepare_data {
    output:
        val geneAnnotation      , emit: geneAnnotation

        path "drop_demo_data-main/Data/rna_bam/*bam"
        file "drop_demo_data-main/Data/*gtf"
    script:
    """
        curl -L https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip -o data.zip && \
            unzip data
    """
}

workflow test_aberrantexpression_counting {
    // get all required input data from demo into corresponding channels  
    (ch_bam, ch_gtf) = prepare_data()
    
    // run the preprocessing only on the gtf channel
    PREPROCESSGENEANNOTATION(ch_gtf)
    
    // run the counting per BAM file
    COUNTREADS(PREPROCESSGENEANNOTATION.out.count_reads, ch_bam.flatten())

    // merge counts
    MERGECOUNTS(COUNTREADS.out.counts.collect(), PREPROCESSGENEANNOTATION.out.count_ranges)

    // filter counts
    FILTERCOUNT(MERGECOUNTS.out.total_counts, PREPROCESSGENEANNOTATION.out.txtdb_out)
}
