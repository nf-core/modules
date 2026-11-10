include { KALLISTOBUSTOOLS_COUNT } from '../../../modules/nf-core/kallistobustools/count/main'
include { KALLISTOBUSTOOLS_REF   } from '../../../modules/nf-core/kallistobustools/ref/main'

workflow FASTQ_QUANT_KALLISTOBUSTOOLS {
    take:
    genome_fasta    // channel: [ meta, fasta ]
    gtf             // channel: [ meta, gtf ]
    kallisto_index  // value channel: prebuilt index
    txp2gene        // path: transcript-to-gene map
    t1c             // path: cDNA capture list
    t2c             // path: intron capture list
    protocol        // val: sequencing technology
    kb_workflow     // val: standard, nac, or lamanno
    ch_fastq        // channel: [ meta, [ reads ] ]

    main:
    assert (txp2gene && kallisto_index) || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

    if (!(txp2gene && kallisto_index)) {
        KALLISTOBUSTOOLS_REF( genome_fasta, gtf, kb_workflow )
        txp2gene = KALLISTOBUSTOOLS_REF.out.t2g.collect()
        kallisto_index = KALLISTOBUSTOOLS_REF.out.index.collect()
        t1c = KALLISTOBUSTOOLS_REF.out.cdna_t2c.ifEmpty{ [] }
        t2c = KALLISTOBUSTOOLS_REF.out.intron_t2c.ifEmpty{ [] }
    } else {
        kallisto_index = kallisto_index.map { index -> [[id: 'kallisto_index'], index] }
    }

    KALLISTOBUSTOOLS_COUNT(
        ch_fastq,
        kallisto_index,
        txp2gene,
        t1c,
        t2c,
        protocol,
        kb_workflow
    )

    ch_raw_counts = KALLISTOBUSTOOLS_COUNT.out.count.map{ meta, kb_dir ->
        if (file("${kb_dir.toUriString()}/counts_unfiltered").exists()) {
            [meta + [input_type: 'raw'], file("${kb_dir.toUriString()}/counts_unfiltered")]
        }
    }
    ch_filtered_counts = KALLISTOBUSTOOLS_COUNT.out.count.map{ meta, kb_dir ->
        if (file("${kb_dir.toUriString()}/counts_filtered").exists()) {
            [meta + [input_type: 'filtered'], file("${kb_dir.toUriString()}/counts_filtered")]
        }
    }

    emit:
    counts          = KALLISTOBUSTOOLS_COUNT.out.count  // channel: [ meta, count directory ]
    counts_raw      = ch_raw_counts                     // channel: [ meta, unfiltered matrix directory ]
    counts_filtered = ch_filtered_counts                // channel: [ meta, filtered matrix directory ]
    txp2gene        = txp2gene                          // path: transcript-to-gene map

}
