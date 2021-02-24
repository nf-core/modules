# Modules Test Data

# Guidelines
TODO: add some guidelines on how to use/ when to add data

Whenever possible, data in the `genomics` directory should be re-used. The `generic` directory is for special cases only. It's use should be minimal.

TODO: remove unecessary files from `generic` as much as possible. Especially:
* fasta
* fastq
* gtf

Those are mainly used for STAR and methyldackel and take up the majority of the space

# Data Description

### genomics

* sarscov2
    * bam:
        * 'sarscov2_paired_aln.bam': sarscvoc2 sequencing reads aligned against GCA_011545545.1_ASM1154554v1_genomic.fasta using minimap2
        * 'test-sc2-artic-v3-sorted-trimmed.bam': sarscov2 reads aligned against MN908947.3 genome
    * bed
        * 'sarscov2.bed': exemplary bed file for MT192765.1 genome
        * 'test-sc2-artic-v3.bed': exemplary bed file for MN908947.3 genome
    * fasta
        * 'GCA_011545545.1_ASM1154554v1_genomic.fasta': MT192765.1 genomem including GATK .dict file
        * 'MN908947.3.fa': MN908947.3 genome
    * fastq
        * 'sarscov2_{1,2}.fastq.gz' sarscov2 paired-end sequencing reads
        * 'sarscov2_b_{1,2}.fastq.gz‘: copies of the above reads
    * gtf
        * 'GCA_011545545.1_ASM1154554v1_genomic.gtf': GTF for MT192765.1 genome
        * 'MN908947.3.gff3': GFF for MN908947.3 genome
        * 'MN908947.3.gff3.gz': bgzipped-version
    * paf
        * 'GCA_011545545.1_ASM1154554v1_cds_from_genomic.paf': PAF file for MT192765.1  genome

### generic
* 'a.gff3.gz': bgzipped gff3 file currently necessary for TABIX test
* bam
    * 'test.paired_end_methylated.sorted.bam': methylated bam file currently necessary for methyldackel
* bedgraph: bedgraph files for seacr
* 'dummy_file.txt': a dummy file for whenever that is required
* fasta: additional fasta files currently necessary for methyldackel and STAR
* fastq: additional fastq files currently necessary for STAR
* gtf: additional gtf file for STAR
* vcf: several VCF files for tools using those, will be removed in the future
* 'test.txt.gar.gz' exemplary tar file for the untar module



