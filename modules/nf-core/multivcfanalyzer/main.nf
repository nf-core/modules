process MULTIVCFANALYZER {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/0713bc8b2176da50f1967ae24dd955f913b7288b966408220596e91c5c14aeeb/data' :
        'community.wave.seqera.io/library/htslib_multivcfanalyzer:6387733614faccad' }"

    input:
    tuple val(meta), path(vcfs)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(snpeff_results)
    tuple val(meta4), path(gff)
    val allele_freqs
    val genotype_quality
    val coverage
    val homozygous_freq
    val heterozygous_freq
    tuple val(meta5), path(gff_exclude)


    output:
    tuple val(meta), path('*fullAlignment.fasta.gz')                                           , emit: full_alignment
    tuple val(meta), path('*info.txt')                                                         , emit: info_txt
    tuple val(meta), path('*snpAlignment.fasta.gz')                                            , emit: snp_alignment
    tuple val(meta), path('*snpAlignmentIncludingRefGenome.fasta.gz')                          , emit: snp_genome_alignment
    tuple val(meta), path('*snpStatistics.tsv')                                                , emit: snpstatistics
    tuple val(meta), path('*snpTable.tsv')                                                     , emit: snptable
    tuple val(meta), path('*snpTableForSnpEff.tsv')                                            , emit: snptable_snpeff
    tuple val(meta), path('*snpTableWithUncertaintyCalls.tsv')                                 , emit: snptable_uncertainty
    tuple val(meta), path('*structureGenotypes.tsv')                                           , emit: structure_genotypes
    tuple val(meta), path('*structureGenotypes_noMissingData-Columns.tsv')                     , emit: structure_genotypes_nomissing
    tuple val(meta), path('*MultiVCFAnalyzer.json')                                            , emit: json
    tuple val("${task.process}"), val('multivcfanalyzer'), eval('multivcfanalyzer -h | head -n 1 | cut -f 3 -d " "') , emit: versions_multivcfanalyzer, topic: versions
    tuple val("${task.process}"), val('tabix'),            eval('tabix -h 2>&1 | grep Version | cut -f 2 -d " "')    , emit: versions_tabix           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: '' // MultiVCFAnalyzer has strict and input ordering and all are mandatory. Deactivating $args to prevent breakage of input
    def args2              = task.ext.args2  ?: ''
    def prefix             = task.ext.prefix ?: "${meta2.id}"
    def cmd_snpeff_results = snpeff_results ? "${snpeff_results}" : "NA"
    def cmd_gff            = gff ? "${gff}" : "NA"
    def cmd_allele_freqs   = allele_freqs ? "T" : "F"
    def cmd_gff_exclude    = gff_exclude ? "${gff}" : "NA"

    """
    multivcfanalyzer \\
        ${cmd_snpeff_results} \\
        ${fasta} \\
        ${cmd_gff} \\
        . \
        ${cmd_allele_freqs}  \\
        ${genotype_quality}  \\
        ${coverage}  \\
        ${homozygous_freq}  \\
        ${heterozygous_freq}  \\
        ${cmd_gff_exclude}  \\
        ${vcfs.sort().join(" ")}

    bgzip \\
        ${args2} \\
        fullAlignment.fasta snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta

    for fn in fullAlignment.fasta.gz info.txt snpAlignment.fasta.gz snpAlignmentIncludingRefGenome.fasta.gz snpStatistics.tsv snpTable.tsv snpTableForSnpEff.tsv snpTableWithUncertaintyCalls.tsv structureGenotypes.tsv structureGenotypes_noMissingData-Columns.tsv MultiVCFAnalyzer.json; do
        mv \${fn} ${prefix}_\${fn}
    done
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    """
    echo "" | gzip > ${prefix}_fullAlignment.fasta.gz
    touch ${prefix}_info.txt
    echo "" | gzip > ${prefix}_snpAlignment.fasta.gz
    echo "" | gzip > ${prefix}_snpAlignmentIncludingRefGenome.fasta.gz
    touch ${prefix}_snpStatistics.tsv
    touch ${prefix}_snpTable.tsv
    touch ${prefix}_snpTableForSnpEff.tsv
    touch ${prefix}_snpTableWithUncertaintyCalls.tsv
    touch ${prefix}_structureGenotypes.tsv
    touch ${prefix}_structureGenotypes_noMissingData-Columns.tsv
    touch ${prefix}_MultiVCFAnalyzer.json
    """
}
