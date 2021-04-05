# Test Data

## Data Access
1. The raw data was retrieved from [this](https://www.ncbi.nlm.nih.gov/bioproject/?term=prjeb39899) project. The two used datasets are [ERR4467723](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=ERR4467723) (tumor) and [ERR4467726](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=ERR4467726) (normal)
2. The data was downloaded using the SRA Toolkit with:
```
prefetch -v <Acc number>
sam-dump <Acc number> | samtools view -bS - > <Acc number>.bam
```
3. Reads mapping to chr 22 were extracted and converted to `fq.gz` using [qbic-pipelines/bamtofastq](https://github.com/qbic-pipelines/bamtofastq)

## Determine region covered by reads

1. Visual inspection
```
chr22   16570000        16610000
```
2. Save in `region.txt`

## VCF reference files:

Following 'reference' vcf files are used. All found in igenomes at `s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/`:
 - dbsnp_146.hg38
 - Mills_and_1000G_gold_standard.indels.hg38
 - gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only

1. Downsize all vcf files to only cover chromosome 22 with:
```
tabix -l dbsnp_146.hg38.vcf.gz | parallel -j 5 'tabix -h dbsnp_146.hg38.vcf.gz {} > {}.vcf'
bgzip dbsnp_chr22.vcf
tabix chr22.vcf
bcftools filter dbsnp_146.hg38.chr22.vcf.gz -r chr22:16570000-16610000 > region_22/dbsnp_146.hg38.chr22_region.vcf
bgzip dbsnp_146.hg38.chr22_region.vcf
tabix dbsnp_146.hg38.chr22_region.vcf.gz
```

2. Manipulate mills & gnomAD file, by changing chr length for chr22 to 40001


## Mapping files:

As base reference `s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes/chr22.fasta` was used.

```
samtools faidx chr22.fasta chr22:16570000-16610000  > region_22/chr22_region.fasta
```

## Sarek pipeline alteration to generate all output files

1. Used release 2.7 container:
2. Add `publishDir` to all UMI related steps
3. Add to mapping process:
```
gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true > ${inputFile1}.fq.gz
```
and `publish` the reads

4. Add `publishDir` to HaplotypeCaller process to publish `.g.vcf` files
5. Run sarek with the following command:
```
nextflow run  ~/.nextflow/assets/nf-core/sarek/main.nf -profile cfc -c sarek.config \
--input 'testdata_dsl2_chr22.tsv' \
--outdir 'results_sarek_22' \
--intervals false  \
--bwa false \
--aligner 'bwa-mem2' \
--igenomes_ignore  \
--save_reference \
--fasta './supportfiles/region_22/chr22_region.fasta' \
--save_bam_mapped \
--genome custom \
--dict false \
--dbsnp './supportfiles/region_22/dbsnp_146.hg38.chr22_region.vcf.gz' \
--dbsnp_index './supportfiles/region_22/dbsnp_146.hg38.chr22_region.vcf.gz.tbi' \
--known_indels './supportfiles/region_22/Mills_and_1000G_gold_standard.indels.hg38.chr22_region.vcf.gz' \
--known_indels_index './supportfiles/region_22/Mills_and_1000G_gold_standard.indels.hg38.chr22_region.vcf.gz.tbi' \
--germline_resource './supportfiles/region_22/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.chr22_region.vcf.gz' \
--germline_resource_index './supportfiles/region_22/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.chr22_region.vcf.gz.tbi' \
--tools 'freebayes,mpileup,msisensor,cnvkit,strelka,HaplotypeCaller,Manta,tiddit' \
--umi --read_structure1 "7M1S+T" --read_structure2 "7M1S+T" \
--max_memory 59.GB \
--max_cpus 19 \
-resume
```

with the following TSV:
```
test    XY      0       testN   1       /sfs/7/workspace/ws/iizha01-dsl2_testdata_human-0/results_chr22/reads/test_normal.1.fq.gz       /sfs/7/workspace/ws/iizha01-dsl2_testdata_human-0/results_chr22/reads/test_normal.2.fq.gz
test    XY      1       testT   2       /sfs/7/workspace/ws/iizha01-dsl2_testdata_human-0/results_chr22/reads/test_tumor.1.fq.gz        /sfs/7/workspace/ws/iizha01-dsl2_testdata_human-0/results_chr22/reads/test_tumor.2.fq.gz
```

## GTF/GFF:

Downloaded the gtf and gff3 files from Ensembl:

1. Download
```
wget http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz
wget http://ftp.ensembl.org/pub/release-103/gff3/homo_sapiens/Homo_sapiens.GRCh38.103.chromosome.22.gff3.gz
```
2. Unzip both with
```
gzip -d
```
3. Copy `region.txt`, change chromosome name to `22`
```
bedtools intersect -a region.txt -b Homo_sapiens.GRCh38.103.chr.gtf -wa -wb > interesting_genes.gtf
```

change 22 to chr22???




ASCAT;
awk '/^chr22/' ../1000G_phase3_GRCh38_maf0.3.loci > 1000G_phase3_GRCh38_maf0.3.chr22.loci

 awk 'NR%2!=0 && NR%3!=0 && NR%4!=0' 1000G_phase3_GRCh38_maf0.3.chr22.loci > 1000G_phase3_GRCh38_maf0.3.chr22_subsample.loci

 awk '$2 == "22"' ../1000G_phase3_GRCh38_maf0.3.loci.gc > ../1000G_phase3_GRCh38_maf0.3.chr22.loci.gc

  awk 'NR%2!=0 && NR%3!=0 && NR%4!=0 && NR%5!=0 && NR%6!=0 && NR%7!=0 && NR%8!=0 && NR%9!=0 && NR%10!=0 && NR%11!=0 && NR%12!=0 && NR%13!=0 && NR%14!=0 && NR%15!=0 && NR%16!=0 && NR%17!=0 && NR%18!=0 && NR%19!=0 && NR%20!=0' ../1000G_phase3_GRCh38_maf0.3.chr22.loci.gc > 1000G_phase3_GRCh38_maf0.3.chr22_subsample.loci.gc


 awk 'NR%2!=0 && NR%3!=0 && NR%4!=0 && NR%5!=0 && NR%6!=0 && NR%7!=0 && NR%8!=0 && NR%9!=0 && NR%10!=0 && NR%11!=0 && NR%12!=0 && NR%13!=0 && NR%14!=0 && NR%15!=0 && NR%16!=0 && NR%17!=0 && NR%18!=0 && NR%19!=0 && NR%20!=0' 1000G_phase3_GRCh38_maf0.3.chr22_subsample.loci.gc > 1000G_phase3_GRCh38_maf0.3.chr22_sub_subsample.loci.gc


## Limitations

Missing files:
1. Contamination tables for Mutect2
2. Reads do not cover chromosome 6
3. Single-end reads
4. Methylated bams
5. Unaligned bams
6. Panel of Normals
7. Ploidy files for ASCAT
8. Mappability files for CONTROLFREEC
9. snpEff & VEP cache
