
include { GALAH                       } from '../../../../modules/nf-core/galah/main.nf'
include { BIOAWK as BIOAWK_CHECKM     } from '../../../../modules/nf-core/bioawk/main.nf'
include { BIOAWK as BIOAWK_GENOMEINFO } from '../../../../modules/nf-core/bioawk/main.nf'
include { GUNZIP                      } from '../../../../modules/nf-core/gunzip/main.nf'

workflow test_galah {

    input = [
        [ id:'test' ], // meta map
        [file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCA_002688505.1_ASM268850v1_genomic.fna.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCF_004296495.1_ASM429649v1_genomic.fna.gz", checkIfExists: true)],
        [], 
        []
    ]

    GALAH ( input )

}

workflow test_galah_genomeinfo {

    genomeinfo = [
        [ id: 'genomeinfo' ],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/magmap/testdata/checkm.lineage_wf.qa_2.tsv", checkIfExists: true)
    ]

    BIOAWK_GENOMEINFO(genomeinfo)

    GUNZIP(BIOAWK_GENOMEINFO.out.output)
    
    ch_genomeinfo = GUNZIP.out.gunzip
        .map { meta, tsv -> [tsv] }

    input = [
        [ id:'test' ], // meta map
        [file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCA_002688505.1_ASM268850v1_genomic.fna.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCF_004296495.1_ASM429649v1_genomic.fna.gz", checkIfExists: true)],
        ch_genomeinfo,
        "genome_info"
    ]



    GALAH ( input )

}
