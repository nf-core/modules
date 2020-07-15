### E. coli paired-end test dataset for Bisulfite-seq applications

The E. coli data set was generated using [Sherman](https://github.com/FelixKrueger/Sherman) as 10,000 reads of paired-end data, with average methylation levels of 80% in CpG context, and 10% in non-CG context. The files can be found in the folder: genaral/fastq/dna, and are called:

```
Ecoli_10K_methylated_R1.fastq.gz
Ecoli_10K_methylated_R2.fastq.gz
```

```bash
Sherman --non_dir --genome /bi/scratch/Genomes/E_coli/ --paired -n 10000 -l 100 --CG 20 --CH 90
```

The data is non-directional, so it should produce roughly 25% mapping to each of the `OT`, `CTOT`, `CTOB` and `OB` strands. Thus, the data can be used for bisulfite mapping in standard (= directional), `--pbat` and `--non_directional` mode.

A test alignment should look roughly like this:

`bismark --genome /bi/scratch/Genomes/E_coli/ -1 Ecoli_10K_methylated_R1.fastq.gz -2 Ecoli_10K_methylated_R2.fastq.gz --non_dir`

```
Bismark report for: Ecoli_10K_methylated_R1.fastq.gz and Ecoli_10K_methylated_R2.fastq.gz (version: v0.22.3)
Bismark was run with Bowtie 2 against the bisulfite genome of /bi/scratch/Genomes/E_coli/ with the specified options: -q --score-min L,0,-0.2 --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500
Option '--non_directional' specified: alignments to all strands were being performed (OT, OB, CTOT, CTOB)

Final Alignment report
======================
Sequence pairs analysed in total:       10000
Number of paired-end alignments with a unique best hit: 9320
Mapping efficiency:     93.2%
Sequence pairs with no alignments under any condition:  0
Sequence pairs did not map uniquely:    680
Sequence pairs which were discarded because genomic sequence could not be extracted:    0

Number of sequence pairs with unique best (first) alignment came from the bowtie output:
CT/GA/CT:       2341    ((converted) top strand)
GA/CT/CT:       2329    (complementary to (converted) top strand)
GA/CT/GA:       2356    (complementary to (converted) bottom strand)
CT/GA/GA:       2294    ((converted) bottom strand)

Final Cytosine Methylation Report
=================================
Total number of C's analysed:   471997

Total methylated C's in CpG context:    111011
Total methylated C's in CHG context:    11923
Total methylated C's in CHH context:    21433
Total methylated C's in Unknown context:        0

Total unmethylated C's in CpG context:  27790
Total unmethylated C's in CHG context:  105877
Total unmethylated C's in CHH context:  193963
Total unmethylated C's in Unknown context:      0

C methylated in CpG context:    80.0%
C methylated in CHG context:    10.1%
C methylated in CHH context:    10.0%
Can't determine percentage of methylated Cs in unknown context (CN or CHN) if value was 0


Bismark completed in 0d 0h 0m 12s
```
