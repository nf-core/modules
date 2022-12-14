#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESSGENEANNOTATION  } from '../../../../../modules/nf-core/aberrantexpression/counting/preprocessgeneannotation'
include { COUNTREADS                } from '../../../../../modules/nf-core/aberrantexpression/counting/countreads'

process prepare_data {
    output:
        path("drop_demo*"), emit: data
    
    script:
    """
        curl -L https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip -o data.zip && \
            unzip data
    """
}

process prepare_data_countreads {
    input:
        val(bam)
        path("count_ranges.rds")

    output:
        val(data)   , emit: data
    
    exec:
        Channel.fromPath(file(bam + ".bam"))
            .set {data}
}

workflow test_aberrantexpression_counting {  
    prepare_data ()
        .map {data  -> [
             "preprocess": [[id: "test"], "$data/Data/gencode_annotation_trunc.gtf"],
             "readcount": "$data/Data/rna_bam/*"]
        }
        .set {input}

    PREPROCESSGENEANNOTATION (
        input.preprocess
    )
    
    prepare_data_countreads (input.readcount, PREPROCESSGENEANNOTATION.out.count_ranges) 
        | COUNTREADS
}