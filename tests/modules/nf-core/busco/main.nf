#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO } from '../../../modules/busco/main.nf'

workflow test_busco_genome_single_fasta {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        ['bacteria_odb10', 'bacteroidetes_odb10'],  // Launch with 'auto' to use --auto-lineage, and specified lineages // 'auto' removed from test due to memory issues
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmpyz_hi62i/busco/
    ├── short_summary.specific.bacteria_odb10.genome.fna.json -> /tmp/tmpza_0dth3/33/7d8c9b2c8931d9ad6a67aa843895e7/short_summary.specific.bacteria_odb10.genome.fna.json
    ├── short_summary.specific.bacteria_odb10.genome.fna.txt -> /tmp/tmpza_0dth3/33/7d8c9b2c8931d9ad6a67aa843895e7/short_summary.specific.bacteria_odb10.genome.fna.txt
    ├── short_summary.specific.bacteroidetes_odb10.genome.fna.json -> /tmp/tmpza_0dth3/6a/e95a0cd21785ce33d63b8f73a68a51/short_summary.specific.bacteroidetes_odb10.genome.fna.json
    ├── short_summary.specific.bacteroidetes_odb10.genome.fna.txt -> /tmp/tmpza_0dth3/6a/e95a0cd21785ce33d63b8f73a68a51/short_summary.specific.bacteroidetes_odb10.genome.fna.txt
    ├── test-bacteria_odb10-busco -> /tmp/tmpza_0dth3/33/7d8c9b2c8931d9ad6a67aa843895e7/test-bacteria_odb10-busco/
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   └── run_bacteria_odb10/
    │   │       ├── busco_sequences/
    │   │       ├── full_table.tsv
    │   │       ├── hmmer_output/
    │   │       ├── missing_busco_list.tsv
    │   │       ├── short_summary.json
    │   │       └── short_summary.txt
    │   └── logs/
    │       └── busco.log
    ├── test-bacteria_odb10-busco.batch_summary.txt -> /tmp/tmpza_0dth3/33/7d8c9b2c8931d9ad6a67aa843895e7/test-bacteria_odb10-busco.batch_summary.txt
    ├── test-bacteroidetes_odb10-busco -> /tmp/tmpza_0dth3/6a/e95a0cd21785ce33d63b8f73a68a51/test-bacteroidetes_odb10-busco/
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   └── run_bacteroidetes_odb10/
    │   │       ├── busco_sequences/
    │   │       ├── full_table.tsv
    │   │       ├── hmmer_output/
    │   │       ├── missing_busco_list.tsv
    │   │       ├── short_summary.json
    │   │       └── short_summary.txt
    │   └── logs/
    │       └── busco.log
    ├── test-bacteroidetes_odb10-busco.batch_summary.txt -> /tmp/tmpza_0dth3/6a/e95a0cd21785ce33d63b8f73a68a51/test-bacteroidetes_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmpza_0dth3/6a/e95a0cd21785ce33d63b8f73a68a51/versions.yml
    */

}

workflow test_busco_genome_multi_fasta {

    input = [
        [ id:'test' ], // meta map
        [
            file( params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true),
            file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['genome_fasta'], checkIfExists: true)
        ]
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmpk19byek7/busco/
    ├── short_summary.specific.bacteria_odb10.genome.fasta.json -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/short_summary.specific.bacteria_odb10.genome.fasta.json
    ├── short_summary.specific.bacteria_odb10.genome.fasta.txt -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/short_summary.specific.bacteria_odb10.genome.fasta.txt
    ├── short_summary.specific.bacteria_odb10.genome.fna.json -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/short_summary.specific.bacteria_odb10.genome.fna.json
    ├── short_summary.specific.bacteria_odb10.genome.fna.txt -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/short_summary.specific.bacteria_odb10.genome.fna.txt
    ├── test-bacteria_odb10-busco -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/test-bacteria_odb10-busco/
    │   ├── genome.fasta/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   └── run_bacteria_odb10/
    │   │       ├── busco_sequences/
    │   │       ├── full_table.tsv
    │   │       ├── hmmer_output/
    │   │       ├── missing_busco_list.tsv
    │   │       ├── short_summary.json
    │   │       └── short_summary.txt
    │   ├── genome.fna/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── prodigal_err.log
    │   │   │   └── prodigal_out.log
    │   │   ├── prodigal_output/
    │   │   │   └── predicted_genes/
    │   │   └── run_bacteria_odb10/
    │   │       ├── busco_sequences/
    │   │       ├── full_table.tsv
    │   │       ├── hmmer_output/
    │   │       ├── missing_busco_list.tsv
    │   │       ├── short_summary.json
    │   │       └── short_summary.txt
    │   └── logs/
    │       └── busco.log
    ├── test-bacteria_odb10-busco.batch_summary.txt -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/test-bacteria_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmplt9fv3tl/15/ff310a16d9ce7ad24e207a05ce718e/versions.yml
    */

}

workflow test_busco_eukaryote_metaeuk {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmpeq4dsir5/busco/
    ├── short_summary.specific.eukaryota_odb10.genome.fasta.json -> /tmp/tmp60hby2pk/6f/529873d91cda6bae3a4a6a21746aee/short_summary.specific.eukaryota_odb10.genome.fasta.json
    ├── short_summary.specific.eukaryota_odb10.genome.fasta.txt -> /tmp/tmp60hby2pk/6f/529873d91cda6bae3a4a6a21746aee/short_summary.specific.eukaryota_odb10.genome.fasta.txt
    ├── test-eukaryota_odb10-busco -> /tmp/tmp60hby2pk/6f/529873d91cda6bae3a4a6a21746aee/test-eukaryota_odb10-busco/
    │   ├── genome.fasta/
    │   │   ├── logs/
    │   │   │   ├── hmmsearch_err.log
    │   │   │   ├── hmmsearch_out.log
    │   │   │   ├── metaeuk_err.log
    │   │   │   └── metaeuk_out.log
    │   │   └── run_eukaryota_odb10/
    │   │       ├── busco_sequences/
    │   │       ├── full_table.tsv
    │   │       ├── hmmer_output/
    │   │       ├── metaeuk_output/
    │   │       ├── missing_busco_list.tsv
    │   │       ├── short_summary.json
    │   │       └── short_summary.txt
    │   └── logs/
    │       └── busco.log
    ├── test-eukaryota_odb10-busco.batch_summary.txt -> /tmp/tmp60hby2pk/6f/529873d91cda6bae3a4a6a21746aee/test-eukaryota_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmp60hby2pk/6f/529873d91cda6bae3a4a6a21746aee/versions.yml
    */

}

workflow test_busco_eukaryote_augustus {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'eukaryota_odb10',
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmp2xqaygjj/busco/
    ├── test-eukaryota_odb10-busco -> /tmp/tmpjqs61x9o/3f/67cc14e873c0ceb45e2a27594d624c/test-eukaryota_odb10-busco/
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
    ├── test-eukaryota_odb10-busco.batch_summary.txt -> /tmp/tmpjqs61x9o/3f/67cc14e873c0ceb45e2a27594d624c/test-eukaryota_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmpjqs61x9o/3f/67cc14e873c0ceb45e2a27594d624c/versions.yml
    */

}

workflow test_busco_protein {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['candidatus_portiera_aleyrodidarum']['genome']['proteome_fasta'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmpzwd5dn56/busco/
    ├── short_summary.specific.bacteria_odb10.proteome.fasta.json -> /tmp/tmpk1nlgbf_/ae/0db07b5cd08fb23d0aba5f134ebbe2/short_summary.specific.bacteria_odb10.proteome.fasta.json
    ├── short_summary.specific.bacteria_odb10.proteome.fasta.txt -> /tmp/tmpk1nlgbf_/ae/0db07b5cd08fb23d0aba5f134ebbe2/short_summary.specific.bacteria_odb10.proteome.fasta.txt
    ├── test-bacteria_odb10-busco -> /tmp/tmpk1nlgbf_/ae/0db07b5cd08fb23d0aba5f134ebbe2/test-bacteria_odb10-busco/
    │   ├── logs/
    │   │   └── busco.log
    │   └── proteome.fasta/
    │       ├── logs/
    │       │   ├── hmmsearch_err.log
    │       │   └── hmmsearch_out.log
    │       └── run_bacteria_odb10/
    │           ├── busco_sequences/
    │           ├── full_table.tsv
    │           ├── hmmer_output/
    │           ├── missing_busco_list.tsv
    │           ├── short_summary.json
    │           └── short_summary.txt
    ├── test-bacteria_odb10-busco.batch_summary.txt -> /tmp/tmpk1nlgbf_/ae/0db07b5cd08fb23d0aba5f134ebbe2/test-bacteria_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmpk1nlgbf_/ae/0db07b5cd08fb23d0aba5f134ebbe2/versions.yml
    */
}
workflow test_busco_transcriptome {

    input = [
        [ id:'test' ], // meta map
        file( params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)
    ]

    BUSCO (
        input,
        'bacteria_odb10',
        [], // Download busco lineage
        [] // No config
    )

    /* Output tree:
    /tmp/tmpitjyvo9g/busco/
    ├── short_summary.specific.bacteria_odb10.test1.contigs.fa.json -> /tmp/tmp6wqi0eyx/4f/ed0b23f0fc807bb68091298845c135/short_summary.specific.bacteria_odb10.test1.contigs.fa.json
    ├── short_summary.specific.bacteria_odb10.test1.contigs.fa.txt -> /tmp/tmp6wqi0eyx/4f/ed0b23f0fc807bb68091298845c135/short_summary.specific.bacteria_odb10.test1.contigs.fa.txt
    ├── test-bacteria_odb10-busco -> /tmp/tmp6wqi0eyx/4f/ed0b23f0fc807bb68091298845c135/test-bacteria_odb10-busco/
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
    ├── test-bacteria_odb10-busco.batch_summary.txt -> /tmp/tmp6wqi0eyx/4f/ed0b23f0fc807bb68091298845c135/test-bacteria_odb10-busco.batch_summary.txt
    └── versions.yml -> /tmp/tmp6wqi0eyx/4f/ed0b23f0fc807bb68091298845c135/versions.yml
    */

}
