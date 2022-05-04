#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf'

// This tests genome decompression, empty input channels and data download
workflow test_busco_genome_single_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        ['bacteria_odb10', 'bacteroidetes_odb10'],
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmpimsfk4sj/busco/
    ├── test-bacteria_odb10-busco -> /tmp/tmp1sz7013h/b7/fdeaab567e1c5bccc475a4c19b8582/test-bacteria_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   ├── run_bacteria_odb10/
    │   │   │   ├── busco_sequences/
    │   │   │   ├── full_table.tsv
    │   │   │   ├── hmmer_output/
    │   │   │   ├── missing_busco_list.tsv
    │   │   │   ├── short_summary.json
    │   │   │   └── short_summary.txt
    │   │   ├── short_summary.specific.bacteria_odb10.genome.fna.json
    │   │   └── short_summary.specific.bacteria_odb10.genome.fna.txt
    │   └── logs/
    │       └── busco.log
    ├── test-bacteroidetes_odb10-busco -> /tmp/tmp1sz7013h/75/0da56f59ee44bd2b85e0172906de49/test-bacteroidetes_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   ├── run_bacteroidetes_odb10/
    │   │   │   ├── busco_sequences/
    │   │   │   ├── full_table.tsv
    │   │   │   ├── hmmer_output/
    │   │   │   ├── missing_busco_list.tsv
    │   │   │   ├── short_summary.json
    │   │   │   └── short_summary.txt
    │   │   ├── short_summary.specific.bacteroidetes_odb10.genome.fna.json
    │   │   └── short_summary.specific.bacteroidetes_odb10.genome.fna.txt
    │   └── logs/
    │       └── busco.log
    └── versions.yml -> /tmp/tmp1sz7013h/b7/fdeaab567e1c5bccc475a4c19b8582/versions.yml
    */

}

workflow test_busco_genome_multi_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
            file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmpt22rjxzq/busco/
    ├── test-bacteria_odb10-busco -> /tmp/tmpfxt64xr_/36/425acbe5e9b27ba0bac8861f735494/test-bacteria_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── genome.fasta/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   ├── run_bacteria_odb10/
    │   │   │   ├── busco_sequences/
    │   │   │   ├── full_table.tsv
    │   │   │   ├── hmmer_output/
    │   │   │   ├── missing_busco_list.tsv
    │   │   │   ├── short_summary.json
    │   │   │   └── short_summary.txt
    │   │   ├── short_summary.specific.bacteria_odb10.genome.fasta.json
    │   │   └── short_summary.specific.bacteria_odb10.genome.fasta.txt
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   ├── run_bacteria_odb10/
    │   │   │   ├── busco_sequences/
    │   │   │   ├── full_table.tsv
    │   │   │   ├── hmmer_output/
    │   │   │   ├── missing_busco_list.tsv
    │   │   │   ├── short_summary.json
    │   │   │   └── short_summary.txt
    │   │   ├── short_summary.specific.bacteria_odb10.genome.fna.json
    │   │   └── short_summary.specific.bacteria_odb10.genome.fna.txt
    │   └── logs/
    │       └── busco.log
    └── versions.yml -> /tmp/tmpfxt64xr_/36/425acbe5e9b27ba0bac8861f735494/versions.yml
    */

}

workflow test_busco_eukaryote_metaeuk {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmp22sf7kg9/busco/
    ├── test-eukaryota_odb10-busco -> /tmp/tmpmic8qsk6/d5/d8cb6681c0fcaa6da34b57ec174d59/test-eukaryota_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── genome.fasta/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── metaeuk_err.log
    │   │   │   └── metaeuk_out.log
    │   │   ├── run_eukaryota_odb10/
    │   │   │   ├── busco_sequences/
    │   │   │   ├── full_table.tsv
    │   │   │   ├── hmmer_output/
    │   │   │   ├── metaeuk_output/
    │   │   │   ├── missing_busco_list.tsv
    │   │   │   ├── short_summary.json
    │   │   │   └── short_summary.txt
    │   │   ├── short_summary.specific.eukaryota_odb10.genome.fasta.json
    │   │   └── short_summary.specific.eukaryota_odb10.genome.fasta.txt
    │   └── logs/
    │       └── busco.log
    └── versions.yml -> /tmp/tmpmic8qsk6/d5/d8cb6681c0fcaa6da34b57ec174d59/versions.yml
    */

}

workflow test_busco_eukaryote_augustus {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmpo77wyvb9/busco/
    ├── test-eukaryota_odb10-busco -> /tmp/tmpshljnwcg/25/9891a19cbabda15a5c10fb5e34987f/test-eukaryota_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── genome.fasta/
    │   │   ├── blast_db/
    │   │   │   ├── genome.fasta.ndb
    │   │   │   ├── genome.fasta.nhr
    │   │   │   ├── genome.fasta.nin
    │   │   │   ├── genome.fasta.not
    │   │   │   ├── genome.fasta.nsq
    │   │   │   ├── genome.fasta.ntf
    │   │   │   └── genome.fasta.nto
    │   │   ├── logs/
    │   │   │   ├── makeblastdb_err.log
    │   │   │   ├── makeblastdb_out.log
    │   │   │   ├── tblastn_err.log
    │   │   │   └── tblastn_out.log
    │   │   └── run_eukaryota_odb10/
    │   │       ├── augustus_output/
    │   │       ├── blast_output/
    │   │       ├── busco_sequences/
    │   │       └── hmmer_output/
    │   └── logs/
    │       └── busco.log
    └── versions.yml -> /tmp/tmpshljnwcg/25/9891a19cbabda15a5c10fb5e34987f/versions.yml
    */

}

workflow test_busco_protein {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmplju98s42/busco/
    ├── test-bacteria_odb10-busco -> /tmp/tmp0oru9_61/9c/e992f5eee84806770002e4510f51cb/test-bacteria_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── logs/
    │   │   └── busco.log
    │   └── proteome.fasta/
    │       ├── logs/
    │       │   ├── hmmsearch_err.log
    │       │   └── hmmsearch_out.log
    │       ├── run_bacteria_odb10/
    │       │   ├── busco_sequences/
    │       │   ├── full_table.tsv
    │       │   ├── hmmer_output/
    │       │   ├── missing_busco_list.tsv
    │       │   ├── short_summary.json
    │       │   └── short_summary.txt
    │       ├── short_summary.specific.bacteria_odb10.proteome.fasta.json
    │       └── short_summary.specific.bacteria_odb10.proteome.fasta.txt
    └── versions.yml -> /tmp/tmp0oru9_61/9c/e992f5eee84806770002e4510f51cb/versions.yml
    */
}
workflow test_busco_transcriptome {

    input = [
        [ id:'test', single_end:false ], // meta map
        file( params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [], // No config
    )

    /* Output tree:
    /tmp/tmp5twpr8o9/busco/
    ├── test-bacteria_odb10-busco -> /tmp/tmp_qyjiads/0d/886515d0f06686b2227517398ef8ce/test-bacteria_odb10-busco/
    │   ├── batch_summary.txt
    │   ├── logs/
    │   │   └── busco.log
    │   └── test1.contigs.fa/
    │       ├── blast_db/
    │       │   ├── test1.contigs.fa.ndb
    │       │   ├── test1.contigs.fa.nhr
    │       │   ├── test1.contigs.fa.nin
    │       │   ├── test1.contigs.fa.not
    │       │   ├── test1.contigs.fa.nsq
    │       │   ├── test1.contigs.fa.ntf
    │       │   └── test1.contigs.fa.nto
    │       ├── logs/
    │       │   ├── hmmsearch_err.log
    │       │   ├── hmmsearch_out.log
    │       │   ├── makeblastdb_err.log
    │       │   ├── makeblastdb_out.log
    │       │   ├── tblastn_err.log
    │       │   └── tblastn_out.log
    │       ├── run_bacteria_odb10/
    │       │   ├── blast_output/
    │       │   ├── busco_sequences/
    │       │   ├── full_table.tsv
    │       │   ├── hmmer_output/
    │       │   ├── missing_busco_list.tsv
    │       │   ├── short_summary.json
    │       │   ├── short_summary.txt
    │       │   └── single_copy_proteins.faa
    │       ├── short_summary.specific.bacteria_odb10.test1.contigs.fa.json
    │       ├── short_summary.specific.bacteria_odb10.test1.contigs.fa.txt
    │       └── translated_proteins/
    │           ├── 1024388at2.faa
    │           ├── 1054741at2.faa
    │           ├── 1093223at2.faa
    │           ├── 1151822at2.faa
    │           ├── 143460at2.faa
    │           ├── 1491686at2.faa
    │           ├── 1504821at2.faa
    │           ├── 1574817at2.faa
    │           ├── 1592033at2.faa
    │           ├── 1623045at2.faa
    │           ├── 1661836at2.faa
    │           ├── 1674344at2.faa
    │           ├── 1698718at2.faa
    │           ├── 1990650at2.faa
    │           ├── 223233at2.faa
    │           ├── 402899at2.faa
    │           ├── 505485at2.faa
    │           ├── 665824at2.faa
    │           ├── 776861at2.faa
    │           ├── 874197at2.faa
    │           ├── 932854at2.faa
    │           └── 95696at2.faa
    └── versions.yml -> /tmp/tmp_qyjiads/0d/886515d0f06686b2227517398ef8ce/versions.yml
    */

}
