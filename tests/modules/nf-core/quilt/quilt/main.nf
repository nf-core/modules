#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QUILT_QUILT                   } from '../../../../../modules/nf-core/quilt/quilt/main.nf'
include { QUILT_QUILT as QUILT_OPTIONAL } from '../../../../../modules/nf-core/quilt/quilt/main.nf'


    // input sequencing data (bam)
    bam = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.haplotagged.1.0.bam', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.ont.1.0.bam',         checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.illumina.1.0.bam',    checkIfExists: true)
    ]

    bai = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.haplotagged.1.0.bam.bai', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.ont.1.0.bam.bai',         checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/NA12878.illumina.1.0.bam.bai',    checkIfExists: true)
    ]

    bamlist = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/bamlist.1.0.txt',                checkIfExists: true)
    ]

    bam_bai_bamlist = [ [ id:"test", chr:"chr20" ], bam, bai, bamlist ]

    // input reference data

    reference_haplotype_file = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.hap.gz', checkIfExists: true)
    ]

    reference_legend_file = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.legend.gz', checkIfExists: true)
    ]

    genetic_map_file = [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/CEU-chr20-final.b38.txt.gz', checkIfExists: true)
    ]

    // parameters

    def chr = "chr20"
    def regions_start = 2000001
    def regions_end = 2100000

    // input channel quilt

    ch_input = [ [ id:"test", chr:"chr20" ], bam, bai, bamlist, reference_haplotype_file, reference_legend_file, chr, regions_start, regions_end , genetic_map_file ]

    // (optional) input truth data

    posfile_phasefile = [
        [ id:'test', chr:"chr20" ], // meta map
        [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.posfile.txt',     checkIfExists: true)
        ],
        [
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/quilt/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.phasefile.txt',   checkIfExists: true)
        ]
    ]

    fasta = [[id:'test'], []]



workflow test_quilt {

    QUILT_QUILT ( ch_input, posfile_phasefile, fasta )
}


workflow test_quilt_no_optional_files {

    posfile = []
    phasefile = []
    posfile_phasefile = [[id: null], posfile, phasefile]
    genetic_map_file = []

    ch_input = [ [ id:"test", chr:"chr20" ], bam, bai, bamlist, reference_haplotype_file, reference_legend_file, chr, regions_start, regions_end, genetic_map_file ]


    QUILT_QUILT ( ch_input, posfile_phasefile, fasta )
}

workflow test_quilt_optional_outputs {

    QUILT_OPTIONAL ( ch_input, posfile_phasefile, fasta )
}
