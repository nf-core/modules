Currently just some notes

* sarscov2_b.bed --> necessary for bedtools intersect, slightly modified only
* sarscov2.fna and sarscov2.fasta made additional *.fasta file because bismark requires this extension
* sarscov2 fastq files with '_b' in the name were created as additional files were needed for the 'cat/fastq' modules
* had to keep the 'a.gff3.gz' file for tabix because the sarscov3 gff wouldn't work
