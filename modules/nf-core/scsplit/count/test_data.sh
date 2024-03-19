#From hadge pipeline
# https://github.com/theislab/hadge/blob/main/test_data/download_data.sh
outputdir="work"
mkdir -p $outputdir && cd $outputdir
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1xl9g3IPFNISa1Z2Uqj3nia4tDWtMDz5H' -O jurkat_293t_exons_only.vcf.withAF.vcf.gz
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1jlhEO1Z7YGYVnxv1YO9arjDwGFeZbfkr' -O jurkat_293t_downsampled_n500_full_bam.bam.bai
FILEID="13CV6CjP9VzmwG5MVHbJiVDMVdiIhGdJB"
FILENAME="jurkat_293t_downsampled_n500_full_bam.bam"
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=$FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$FILEID" -O $FILENAME && rm -rf /tmp/cookies.txt
