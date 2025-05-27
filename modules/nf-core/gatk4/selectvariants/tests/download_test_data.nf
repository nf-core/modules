process DOWNLOAD_TEST_DATA {
    output:
    path "*"

    script:
    """
    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gvcf/test.genome.vcf.gz
    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/gvcf/test.genome.vcf.gz.tbi
    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta
    wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta.fai
    """
}

workflow {
    DOWNLOAD_TEST_DATA()
}
