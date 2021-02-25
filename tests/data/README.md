# Modules Test Data

This directory contains all data used for the individual module tests. It is currently organised in `genomics` and `generic`. The former contains all typical data required for genomics modules, such as fasta, fastq and bam files. Every folder in `genomics` corresponds to a single organisms. Any other data is stored in `generic`. This contains files that currently cannot be associated to a genomics category, but also depreciated files which will be removed in the future and exchanged by files in `genomics`.

When adding a new module, please check carefully whether the data necessary for the tests exists already in `tests/data/genomics`. If you can't find the data, please ask about it in the slack #modules channel.

## Data Description

### genomics

* sarscov2
    * bam:
        * 'test_paired_end.bam': sarscvoc2 sequencing reads aligned against test_genomic.fasta using minimap2
        * 'test_paired_end.sorted.bam': sorted version of the above bam file
        * 'test_paired_end.bam.sorted.bam.bai': bam index for the sorted bam file
        * 'test_single_end.bam': alignment (unsorted) of the 'test_1.fastq.gz' reads against test_genomic.fasta using minimap2
    * bed
        * 'test.bed': exemplary bed file for the MT192765.1 genome (fasta/test_genomic.fasta)
        * 'test.2.bed': slightly modified copy of the above file
        * 'test.bed.gz': gzipped version
        * 'test.genome.sizes': genome size for the MT192765.1 genome
    * fasta
        * 'test_genomic.fasta': MT192765.1 genomem including (GCA_011545545.1_ASM1154554v1)
        * 'test_genomic.dict': GATK dict for 'test_genomic.fasta'
        * 'test_genomic.fasta.fai': fasta index for 'test_genomic.fasta'
        * 'test_cds_from_genomic.fasta': coding sequencing from MT192765.1 genome (transcripts)
    * fastq
        * 'test_{1,2}.fastq.gz' sarscov2 paired-end sequencing reads
        * 'test_{1,2}.2.fastq.gz‘: copies of the above reads
    * gtf
        * 'test_genomic.gtf': GTF for MT192765.1 genome
        * 'test_genomic.gff3': GFF for MT192765.1 genome
        * 'test_genomic.gff3.gz': bgzipped-version
    * paf
        * 'test_cds_from_genomic.paf': PAF file for MT192765.1  genome

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



