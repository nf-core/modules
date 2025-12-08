# This dict provides 'asset packages', which specify recipes (commands) to
# build assets. Each package can produce one or more assets, which are encoded
# as relative paths. The package name is often the same as the asset name but
# it need not be.

# These building recipes should make use of arguments that are auto-populated,
# or user-provided. The auto-populated arguments are:
# - {genome}
# - {asset_outfolder}
#   In addition to these, the recipe should refer in the
#   same way, {var}, to any variables required to be provided, which will be
#   provided via the CLI. These should be listed as 'required_inputs' and
#   will be checked for existence before the commands are executed.

DESC = "description"
ASSET_DESC = "asset_description"
ASSETS = "assets"
PTH = "path"
REQ_FILES = "required_files"
REQ_ASSETS = "required_assets"
REQ_PARAMS = "required_parameters"
CONT = "container"
CMD_LST = "command_list"
KEY = "key"
DEFAULT = "default"

RECIPE_CONSTS = [
    "DESC",
    "ASSET_DESC",
    "ASSETS",
    "PTH",
    "REQ_FILES",
    "REQ_ASSETS",
    "CONT",
    "CMD_LST",
    "KEY",
    "DEFAULT",
]

asset_build_packages = {
    "fasta": {
        DESC: "DNA sequences in the FASTA format, indexed FASTA (produced with samtools index) and chromosome sizes file",
        ASSETS: {
            "fasta": "{genome}.fa",
            "fai": "{genome}.fa.fai",
            "chrom_sizes": "{genome}.chrom.sizes",
        },
        REQ_FILES: [{KEY: "fasta", DESC: "gzipped fasta file"}],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {fasta} {asset_outfolder}/{genome}.fa.gz",
            "gzip -df {asset_outfolder}/{genome}.fa.gz",
            "samtools faidx {asset_outfolder}/{genome}.fa",
            "cut -f 1,2 {asset_outfolder}/{genome}.fa.fai > {asset_outfolder}/{genome}.chrom.sizes",
        ],
    },
    "fasta_txome": {
        DESC: "cDNA sequences in the FASTA format, indexed FASTA (produced with samtools index) and chromosome sizes file",
        ASSETS: {
            "fasta_txome": "{genome}.fa",
            "fai": "{genome}.fa.fai",
            "chrom_sizes": "{genome}.chrom.sizes",
        },
        REQ_FILES: [{KEY: "fasta", DESC: "gzipped fasta file"}],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {fasta} {asset_outfolder}/{genome}.fa.gz",
            "gzip -df {asset_outfolder}/{genome}.fa.gz",
            "samtools faidx {asset_outfolder}/{genome}.fa",
            "cut -f 1,2 {asset_outfolder}/{genome}.fa.fai > {asset_outfolder}/{genome}.chrom.sizes",
        ],
    },
    "dbnsfp": {
        DESC: "A database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome (Gencode release 29/Ensembl 94)",
        ASSETS: {
            "dbnsfp": "{genome}_dbNSFP.txt.gz",
            "tabix": "{genome}_dbNSFP.txt.gz.tbi",
        },
        REQ_FILES: [{KEY: "dbnsfp", DESC: "zipped dbNSFP database file"}],
        REQ_ASSETS: [],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            }
        ],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {dbnsfp} {asset_outfolder}/{genome}.zip",
            "unzip {asset_outfolder}/{genome}.zip -d {asset_outfolder}",
            "gunzip -v {asset_outfolder}/*variant.chr*.gz",
            "head -n1 {asset_outfolder}/dbNSFP*_variant.chr1 > {asset_outfolder}/{genome}_dbNSFP.txt",
            "cat {asset_outfolder}/dbNSFP*variant.chr* | grep -v '#' >> {asset_outfolder}/{genome}_dbNSFP.txt",
            "rm {asset_outfolder}/dbNSFP*_variant.chr*",
            "bgzip -@ {threads} {asset_outfolder}/{genome}_dbNSFP.txt",
            "tabix -s 1 -b 2 -e 2 {asset_outfolder}/{genome}_dbNSFP.txt.gz",
            "rm `find {asset_outfolder} -type f -not -path '{asset_outfolder}/_refgenie_build*' -not -path '{asset_outfolder}/{genome}_dbNSFP.txt.*'`",
        ],
    },
    "dbsnp": {
        DESC: "The database of single nucleotide polymorphisms (SNPs) and multiple small-scale variations that include insertions/deletions, microsatellites, and non-polymorphic variants",
        ASSETS: {"dbsnp": "{genome}_dbSNP.gz", "tabix": "{genome}_dbSNP.gz.tbi"},
        REQ_FILES: [
            {KEY: "dbsnp_vcf", DESC: "SNP database file in Variant Call Format (VCF)"},
            {KEY: "dbsnp_tbi", DESC: "tabix index of the dbsnp.vcf file"},
        ],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {dbsnp_vcf} {asset_outfolder}/{genome}_dbSNP.gz",
            "cp {dbsnp_tbi} {asset_outfolder}/{genome}_dbSNP.gz.tbi",
        ],
    },
    "bowtie2_index": {
        DESC: "Genome index for bowtie, produced with bowtie-build",
        ASSETS: {"bowtie2_index": "{genome}"},
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: ["bowtie2-build {fasta} {asset_outfolder}/{genome}"],
    },
    "bwa_index": {
        DESC: "Genome index for Burrows-Wheeler Alignment Tool, produced with bwa index",
        ASSETS: {"bwa_index": "{genome}.fa"},
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bwa index {asset_outfolder}/{genome}.fa",
        ],
    },
    "hisat2_index": {
        DESC: "Genome index for HISAT2, produced with hisat2-build",
        ASSETS: {"hisat2_index": "{genome}"},
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: ["hisat2-build {fasta} {asset_outfolder}/{genome}"],
    },
    "bismark_bt2_index": {
        DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie2",
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {"bismark_bt2_index": "."},
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bismark_genome_preparation --bowtie2 {asset_outfolder}",
        ],
    },
    "bismark_bt1_index": {
        DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie1",
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {"bismark_bt1_index": "."},
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bismark_genome_preparation {asset_outfolder}",
        ],
    },
    "kallisto_index": {
        DESC: "Genome index for kallisto, produced with kallisto index",
        REQ_FILES: [],
        REQ_ASSETS: [
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for transcriptome"}
        ],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {"kallisto_index": "."},
        CMD_LST: [
            "kallisto index -i {asset_outfolder}/{genome}_kallisto_index.idx {fasta}"
        ],
    },
    "salmon_index": {
        DESC: "Transcriptome index for salmon, produced with salmon index",
        REQ_FILES: [],
        REQ_ASSETS: [
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for transcriptome"}
        ],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            },
            {
                KEY: "kmer",
                DEFAULT: "31",
                DESC: "The length of kmer to use to create the indices",
            },
        ],
        CONT: "combinelab/salmon",
        ASSETS: {"salmon_index": "."},
        CMD_LST: [
            "salmon index -t {fasta} -i {asset_outfolder} -k {kmer} -p {threads}"
        ],
    },
    "salmon_sa_index": {
        DESC: "Transcriptome index for salmon, produced with salmon index using selective alignment method. Improves quantification accuracy compared to the regular index.",
        REQ_FILES: [],
        REQ_ASSETS: [
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"},
            {
                KEY: "fasta_txome",
                DEFAULT: "fasta_txome",
                DESC: "fasta asset for transcriptome",
            },
        ],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            },
            {
                KEY: "kmer",
                DEFAULT: "31",
                DESC: "The length of kmer to use to create the indices",
            },
        ],
        CONT: "combinelab/salmon",
        ASSETS: {"salmon_sa_index": "."},
        CMD_LST: [
            "grep '^>' {fasta} | cut -d ' ' -f 1 > {asset_outfolder}/decoys.txt",
            "sed -i.bak -e 's/>//g' {asset_outfolder}/decoys.txt",
            "rm {asset_outfolder}/decoys.txt.bak",
            "cat {fasta_txome} {fasta} > {asset_outfolder}/gentrome.fa",
            "salmon index -t {asset_outfolder}/gentrome.fa -d {asset_outfolder}/decoys.txt -i {asset_outfolder} -k {kmer} -p {threads}",
            "rm {asset_outfolder}/gentrome.fa {asset_outfolder}/decoys.txt",
        ],
    },
    "salmon_partial_sa_index": {
        DESC: "Transcriptome index for salmon, produced with salmon index using partial selective alignment method. Preparation includes transcriptome mapping to the genome and extraction of the relevant portion out from the genome and indexing it along with the transcriptome. Recipe source -- https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh",
        REQ_FILES: [],
        REQ_ASSETS: [
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"},
            {
                KEY: "fasta_txome",
                DEFAULT: "fasta_txome",
                DESC: "fasta asset for transcriptome",
            },
            {
                KEY: "gtf",
                DEFAULT: "ensembl_gtf",
                DESC: "GTF file for exonic features extraction",
            },
        ],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            },
            {
                KEY: "kmer",
                DEFAULT: "31",
                DESC: "The length of kmer to use to create the indices",
            },
        ],
        CONT: "combinelab/salmon",
        ASSETS: {"salmon_partial_sa_index": "."},
        CMD_LST: [
            "gunzip -c {gtf} > {asset_outfolder}/{genome}.gtf",
            "awk -v OFS='\t' '{{if ($3==\"exon\") {{print $1,$4,$5}}}}' {asset_outfolder}/{genome}.gtf > {asset_outfolder}/exons.bed",
            "bedtools maskfasta -fi {fasta} -bed {asset_outfolder}/exons.bed -fo {asset_outfolder}/reference.masked.genome.fa",
            "mashmap -r {asset_outfolder}/reference.masked.genome.fa -q {fasta_txome} -t {threads} --pi 80 -s 500 -o {asset_outfolder}/mashmap.out",
            "awk -v OFS='\t' '{{print $6,$8,$9}}' {asset_outfolder}/mashmap.out | sort -k1,1 -k2,2n - > {asset_outfolder}/genome_found.sorted.bed",
            "bedtools merge -i {asset_outfolder}/genome_found.sorted.bed > {asset_outfolder}/genome_found_merged.bed",
            "bedtools getfasta -fi {asset_outfolder}/reference.masked.genome.fa -bed {asset_outfolder}/genome_found_merged.bed -fo {asset_outfolder}/genome_found.fa",
            'awk \'{{a=$0; getline;split(a, b, ":");  r[b[1]] = r[b[1]]""$0}} END {{ for (k in r) {{ print k"\\n"r[k] }} }}\' {asset_outfolder}/genome_found.fa > {asset_outfolder}/decoy.fa',
            "cat {fasta_txome} {asset_outfolder}/decoy.fa > {asset_outfolder}/gentrome.fa",
            "grep '>' {asset_outfolder}/decoy.fa | awk '{{print substr($1,2); }}' > {asset_outfolder}/decoys.txt",
            "rm {asset_outfolder}/exons.bed {asset_outfolder}/reference.masked.genome.fa {asset_outfolder}/mashmap.out {asset_outfolder}/genome_found.sorted.bed {asset_outfolder}/genome_found_merged.bed {asset_outfolder}/genome_found.fa {asset_outfolder}/decoy.fa {asset_outfolder}/reference.masked.genome.fa.fai",
            "salmon index -t {asset_outfolder}/gentrome.fa -d {asset_outfolder}/decoys.txt -i {asset_outfolder} -k {kmer} -p {threads}",
        ],
    },
    "tgMap": {
        DESC: "Transcript to gene map file, containing two columns mapping of each transcript present in the reference to the corresponding gene.",
        REQ_FILES: [],
        REQ_ASSETS: [
            {
                KEY: "salmon_partial_sa_index",
                DEFAULT: "salmon_partial_sa_index",
                DESC: "partial salmon index asset",
            }
        ],
        REQ_PARAMS: [],
        ASSETS: {"tgMap": "{genome}_txp2gene.tsv"},
        CMD_LST: [
            "grep '^>' {salmon_partial_sa_index}/gentrome.fa | cut -d ' ' -f 1,7 | tr -s ' ' '\\t' | sed 's/[>'gene_symbol:']//g' > {asset_outfolder}/{genome}_txp2gene.tsv",
        ],
    },
    "epilog_index": {
        DESC: "Genome index for CpG sites, produced by the epilog DNA methylation caller",
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [
            {
                KEY: "context",
                DEFAULT: "CG",
                DESC: "Substring to index. One or more space-separated strings to index. e.g. 'CG' or 'CG CA CT CC'",
            }
        ],
        CONT: "databio/refgenie",
        ASSETS: {"epilog_index": "{genome}_{context}.tsv.gz"},
        CMD_LST: [
            "epilog index -- --infile {fasta} --outfile {asset_outfolder}/{genome}_{context}.tsv --contexts {context}",
            "bgzip {asset_outfolder}/{genome}_{context}.tsv",
            "tabix -s 1 -b 2 -e 2 {asset_outfolder}/{genome}_{context}.tsv.gz",
        ],
    },
    "star_index": {
        DESC: "Genome index for STAR RNA-seq aligner, produced with STAR --runMode genomeGenerate",
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            }
        ],
        CONT: "databio/refgenie",
        ASSETS: {"star_index": "."},
        CMD_LST: [
            "mkdir -p {asset_outfolder}",
            "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {asset_outfolder} --genomeFastaFiles {fasta}",
        ],
    },
    "gencode_gtf": {
        DESC: "GTF annotation asset which provides access to all annotated transcripts which make up an Ensembl gene set.",
        REQ_FILES: [
            {
                KEY: "gencode_gtf",
                DESC: "Annotation file in Gene Transfer Format (GTF) from Gencode",
            }
        ],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {"gencode_gtf": "{genome}.gtf.gz"},
        CMD_LST: ["cp {gencode_gtf} {asset_outfolder}/{genome}.gtf.gz"],
    },
    "ensembl_gtf": {
        DESC: "Ensembl GTF, TSS, and gene body annotation",
        REQ_FILES: [
            {
                KEY: "ensembl_gtf",
                DESC: "Annotation file in Gene Transfer Format (GTF) from Ensembl",
            }
        ],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "ensembl_gtf": "{genome}.gtf.gz",
            "ensembl_tss": "{genome}_ensembl_TSS.bed",
            "ensembl_gene_body": "{genome}_ensembl_gene_body.bed",
        },
        CMD_LST: [
            "cp {ensembl_gtf} {asset_outfolder}/{genome}.gtf.gz",
            'gzip -dcf {asset_outfolder}/{genome}.gtf.gz | grep \'exon_number "1";\' | sed \'s/^/chr/\' | awk -v OFS=\'\t\' \'{{print $1, $4, $5, $20, $14, $7}}\' | sed \'s/";//g\' | sed \'s/"//g\' | awk \'{{if($6=="+"){{print $1"\t"$2+20"\t"$2+120"\t"$4"\t"$5"\t"$6}}else{{print $1"\t"$3-120"\t"$3-20"\t"$4"\t"$5"\t"$6}}}}\' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_ensembl_TSS.bed',
            'gzip -dcf {asset_outfolder}/{genome}.gtf.gz | awk \'$3 == "gene"\' | sed \'s/^/chr/\' | awk -v OFS=\'\t\' \'{{print $1,$4,$5,$14,$6,$7}}\' | sed \'s/";//g\' | sed \'s/"//g\' | awk \'$4!="Metazoa_SRP"\' | awk \'$4!="U3"\' | awk \'$4!="7SK"\'  | awk \'($3-$2)>200\' | awk \'{{if($6=="+"){{print $1"\t"$2+500"\t"$3"\t"$4"\t"$5"\t"$6}}else{{print $1"\t"$2"\t"$3-500"\t"$4"\t"$5"\t"$6}}}}\' | awk \'$3>$2\' | LC_COLLATE=C sort -k4 -u > {asset_outfolder}/{genome}_ensembl_gene_body.bed',
        ],
    },
    "ensembl_rb": {
        DESC: "A regulatory annotation file",
        REQ_FILES: [
            {
                KEY: "gff",
                DESC: "Regulatory build annotation file in Gene Feature Format (GFF) from Ensembl",
            }
        ],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {"ensembl_rb": "{genome}.gff.gz"},
        CMD_LST: ["cp {gff} {asset_outfolder}/{genome}.gff.gz"],
    },
    "refgene_anno": {
        DESC: "gene, TSS, exon, intron, and premature mRNA annotation files",
        REQ_FILES: [{KEY: "refgene", DESC: "gzipped RefGene database annotation file"}],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "refgene_anno": "{genome}_refGene.txt.gz",
            "refgene_tss": "{genome}_TSS.bed",
            "refgene_exon": "{genome}_exons.bed",
            "refgene_intron": "{genome}_introns.bed",
            "refgene_pre_mRNA": "{genome}_pre-mRNA.bed",
        },
        CMD_LST: [
            "cp {refgene} {asset_outfolder}/{genome}_refGene.txt.gz",
            'gzip -dcf {asset_outfolder}/{genome}_refGene.txt.gz | awk \'{{if($4=="+"){{print $3"\t"$5"\t"$5"\t"$13"\t.\t"$4}}else{{print $3"\t"$6"\t"$6"\t"$13"\t.\t"$4}}}}\' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_TSS.bed',
            "gzip -dcf {asset_outfolder}/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -v OFS='\t' '{{ n = split($10, a, \",\"); split($11, b, \",\"); for(i=1; i<n; ++i) print $3, a[i], b[i], $13, i, $4 }}' | awk -v OFS='\t' '$6==\"+\" && $5!=1 {{print $0}} $6==\"-\" {{print $0}}' | awk '$4!=prev4 && prev6==\"-\" {{prev4=$4; prev6=$6; delete line[NR-1]; idx-=1}} {{line[++idx]=$0; prev4=$4; prev6=$6}} END {{for (x=1; x<=idx; x++) print line[x]}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dcf {asset_outfolder}/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -F'\t' '{{ exonCount=int($9);split($10,exonStarts,\"[,]\"); split($11,exonEnds,\"[,]\"); for(i=1;i<exonCount;i++) {{printf(\"%s\\t%s\\t%s\\t%s\\t%d\\t%s\\n\",$3,exonEnds[i],exonStarts[i+1],$13,($3==\"+\"?i:exonCount-i),$4);}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            'gzip -dcf {asset_outfolder}/{genome}_refGene.txt.gz | grep \'cmpl\' | awk  \'{{print $3"\t"$5"\t"$6"\t"$13"\t.\t"$4}}\' | LC_COLLATE=C sort -k1,1 -k2,2n -u >  {asset_outfolder}/{genome}_pre-mRNA.bed',
        ],
    },
    "suffixerator_index": {
        DESC: "Enhanced suffix array index for genomes using gt (GenomeTools) suffixerator program",
        REQ_PARAMS: [
            {
                KEY: "memlimit",
                DEFAULT: "8GB",
                DESC: "The maximum amount of memory available to be used during index construction.",
            }
        ],
        REQ_FILES: [],
        REQ_ASSETS: [{KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"}],
        CONT: "databio/refgenie",
        ASSETS: {"esa": "{genome}.sft"},
        CMD_LST: [
            "gt suffixerator -dna -pl -tis -suf -lcp -v -showprogress -memlimit {memlimit} -db {fasta} -indexname {asset_outfolder}/{genome}.sft"
        ],
    },
    "tallymer_index": {
        DESC: "Indexed k-mers for a given enhanced suffix array at a fixed value of k",
        REQ_PARAMS: [
            {KEY: "mersize", DEFAULT: "30", DESC: "The mer size."},
            {
                KEY: "minocc",
                DEFAULT: "2",
                DESC: "The minimum occurrence number for the mers to index.",
            },
        ],
        REQ_FILES: [],
        REQ_ASSETS: [
            {
                KEY: "esa",
                DEFAULT: "suffixerator_index",
                DESC: "enhanced suffix array index for genome",
            },
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"},
        ],
        CONT: "databio/refgenie",
        ASSETS: {
            "tindex": "{genome}.tal_{mersize}",
            "search_file": "{genome}.tal_{mersize}.gtTxt",
        },
        CMD_LST: [
            "gt tallymer mkindex -v -counts -pl -mersize {mersize} -minocc {minocc} -indexname {asset_outfolder}/{genome}.tal_{mersize} -esa {esa}/{genome}.sft",
            "gt tallymer search -output qseqnum qpos -strand fp -tyr {asset_outfolder}/{genome}.tal_{mersize} -q {fasta} > {asset_outfolder}/{genome}.tal_{mersize}.gtTxt",
        ],
    },
    "feat_annotation": {
        DESC: "Combined genomic feature annotation created using an Ensembl GTF annotation asset and an Ensembl regulatory build annotation asset",
        ASSETS: {
            "feat_annotation": "{genome}_annotations.bed.gz",
        },
        REQ_FILES: [],
        REQ_ASSETS: [
            {
                KEY: "ensembl_gtf",
                DEFAULT: "ensembl_gtf",
                DESC: "Annotation file in Gene Transfer Format (GTF) from Ensembl",
            },
            {
                KEY: "ensembl_rb",
                DEFAULT: "ensembl_rb",
                DESC: "Regulatory annotation file in General Feature Format (GTF) from Ensembl",
            },
        ],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "gzip -dcf {ensembl_gtf} | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"Exon\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dcf {ensembl_gtf} | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{ split($20, a, \"\\\"\"); print \"chr\"$1, $4-1, $5, a[2], $6, $7}}' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u | awk 'seen[$4]++ && seen[$4] > 1' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3nr | env LC_COLLATE=C sort -k1,1 -k2,2n -u | env LC_COLLATE=C sort -k1,1 -k3,3n -u | awk -v OFS='\t' '{{if($4==prev4){{new2=prev3+1;}} {{prev4=$4; prev3=$3; print $1, new2, $2-1, \"Intron\", $5, $6}}}}' | awk -F'\t' '$2' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dcf {ensembl_gtf} | awk '$3==\"three_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"3'\\'' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_3utr.bed",
            "gzip -dcf {ensembl_gtf} | awk '$3==\"five_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"5'\\'' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_5utr.bed",
            "gzip -dcf {ensembl_rb} | awk '$3==\"promoter\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter.bed",
            "gzip -dcf {ensembl_rb} | awk '$3==\"promoter_flanking_region\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter Flanking Region\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter_flanking.bed",
            "gzip -dcf {ensembl_rb} | awk '$3==\"enhancer\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Enhancer\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_enhancer.bed",
            "cat {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed | awk -F'\t' '!seen[$1, $2, $3]++' > {asset_outfolder}/{genome}_annotations.bed",
            "rm -f {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed",
            "gzip -f {asset_outfolder}/{genome}_annotations.bed",
        ],
    },
    "cellranger_reference": {
        DESC: "Cell Ranger custom genome reference for read alignment and gene expression quantification",
        ASSETS: {
            "cellranger_reference": "ref",
        },
        REQ_FILES: [],
        REQ_ASSETS: [
            {
                KEY: "gtf",
                DEFAULT: "gencode_gtf",
                DESC: "Annotation file in Gene Transfer Format (GTF) from Gencode",
            },
            {KEY: "fasta", DEFAULT: "fasta", DESC: "fasta asset for genome"},
        ],
        REQ_PARAMS: [
            {
                KEY: "threads",
                DEFAULT: "8",
                DESC: "Number of threads to use for parallel computing",
            }
        ],
        CONT: "databio/refgenie",
        CMD_LST: [
            "gunzip {gtf} -c > {asset_outfolder}/{genome}.gtf",
            "cellranger mkgtf {asset_outfolder}/{genome}.gtf {asset_outfolder}/{genome}_filtered.gtf",
            "rm {asset_outfolder}/{genome}.gtf",
            "cd {asset_outfolder}; cellranger mkref --genome=ref --fasta={fasta} --genes={asset_outfolder}/{genome}_filtered.gtf --nthreads={threads}",
        ],
    },
    "blacklist": {
        DESC: "Atypical, unstructured, or high signal genomic regions present in next-generation sequencing experiments (e.g. from ENCODE)",
        ASSETS: {
            "blacklist": "{genome}_blacklist.bed.gz",
        },
        REQ_FILES: [{KEY: "blacklist", DESC: "gzipped blacklist file"}],
        REQ_ASSETS: [],
        REQ_PARAMS: [],
        CONT: "databio/refgenie",
        CMD_LST: ["cp {blacklist} {asset_outfolder}/{genome}_blacklist.bed.gz"],
    },
}
