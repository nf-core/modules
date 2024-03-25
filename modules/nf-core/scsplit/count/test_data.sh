#From hadge pipeline
# https://github.com/theislab/hadge/blob/main/test_data/download_data.sh
outputdir="work"
mkdir -p $outputdir && cd $outputdir
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1xl9g3IPFNISa1Z2Uqj3nia4tDWtMDz5H' -O jurkat_293t_exons_only.vcf.withAF.vcf.gz
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1jlhEO1Z7YGYVnxv1YO9arjDwGFeZbfkr' -O jurkat_293t_downsampled_n500_full_bam.bam.bai
FILEID="13CV6CjP9VzmwG5MVHbJiVDMVdiIhGdJB"
FILENAME="jurkat_293t_downsampled_n500_full_bam.bam"
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=$FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$FILEID" -O $FILENAME && rm -rf /tmp/cookies.txt
wget --no-check-certificate https://figshare.com/ndownloader/files/41773428 -O barcodes.tsv


gunzip jurkat_293t_exons_only.vcf.withAF.vcf.gz

# Download subset reference genome
wget --no-check-certificate https://figshare.com/ndownloader/files/43102459 -O genome_chr1.fa
wget --no-check-certificate https://figshare.com/ndownloader/files/43102453 -O genome_chr1.fa.fai
# source: http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-3.0.0.tar.gz

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y" -O common_variants_hg19.vcf && rm -rf /tmp/cookies.txt
wget https://master.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz

freebayes -f genome_chr1.fa jurkat_293t_downsampled_n500_full_bam.bam > snv.vcf
