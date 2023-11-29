#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTA_NEWICK_EPANG_GAPPA } from '../../../../subworkflows/nf-core/fasta_newick_epang_gappa/main.nf'

workflow test_fasta_newick_epang_gappa_hmmer {
    
    ch_input = Channel.of(
        [
            meta: [ id: "hmmer" ],
            data: [
                alignmethod:  'hmmer',
                queryseqfile: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_3_sequences.faa", checkIfExists: true),
                refseqfile:   file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_seed.alnfaa", checkIfExists: true),
                refphylogeny: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_seed.ft.LGCAT.newick", checkIfExists: true),
                model:        "LG",
                taxonomy:     file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/gappa/gappa_taxonomy.tsv", checkIfExists: true)
            ]
        ]
    )

    FASTA_NEWICK_EPANG_GAPPA ( ch_input )
}

workflow test_fasta_newick_epang_gappa_mafft {
    
    ch_input = Channel.of(
        [
            meta: [ id: "mafft" ],
            data: [
                alignmethod:  'mafft',
                queryseqfile: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_3_sequences.faa", checkIfExists: true),
                refseqfile:   file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_seed.alnfaa", checkIfExists: true),
                refphylogeny: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/PF14720_seed.ft.LGCAT.newick", checkIfExists: true),
                model:        "LG",
                taxonomy:     file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/gappa/gappa_taxonomy.tsv", checkIfExists: true)
            ]
        ]
    )

    FASTA_NEWICK_EPANG_GAPPA ( ch_input )
}

workflow test_fasta_newick_epang_gappa_nucl_hmmer {
    
    ch_input = Channel.of(
        [
            meta: [ id: "nucl_hmmer" ],
            data: [
                alignmethod:  'hmmer',
                queryseqfile: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/cyn_syn.fna", checkIfExists: true),
                refseqfile:   file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/cyanos_16s.alnfna", checkIfExists: true),
                refphylogeny: file("https://github.com/nf-core/test-datasets/raw/phyloplace/testdata/cyanos_16s.newick", checkIfExists: true),
                model:        "GTR+F+I+I+R3",
                taxonomy:     []
            ]
        ]
    )

    FASTA_NEWICK_EPANG_GAPPA ( ch_input )
}
