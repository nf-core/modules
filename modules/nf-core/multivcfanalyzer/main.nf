process MULTIVCFANALYZER {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83ae9e52b1c7d8dfec31c1e045908153c9c575b212f0d89c56885315f51d98d3/data' :
        'community.wave.seqera.io/library/htslib_multivcfanalyzer:71bddf4ac37aff5b' }"

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
    tuple val(meta), path('fullAlignment.fasta.gz')                                           , emit: full_alignment
    tuple val(meta), path('info.txt')                                                         , emit: info_txt
    tuple val(meta), path('snpAlignment.fasta.gz')                                            , emit: snp_alignment
    tuple val(meta), path('snpAlignmentIncludingRefGenome.fasta.gz')                          , emit: snp_genome_alignment
    tuple val(meta), path('snpStatistics.tsv')                                                , emit: snpstatistics
    tuple val(meta), path('snpTable.tsv')                                                     , emit: snptable
    tuple val(meta), path('snpTableForSnpEff.tsv')                                            , emit: snptable_snpeff
    tuple val(meta), path('snpTableWithUncertaintyCalls.tsv')                                 , emit: snptable_uncertainty
    tuple val(meta), path('structureGenotypes.tsv')                                           , emit: structure_genotypes
    tuple val(meta), path('structureGenotypes_noMissingData-Columns.tsv')                     , emit: structure_genotypes_nomissing
    tuple val(meta), path('MultiVCFAnalyzer.json')                                            , emit: json
    tuple val("${task.process}"), val('multivcfanalyzer'), eval('multivcfanalyzer -h | head -n 1 | cut -f 3 -d " "') , emit: versions_multivcfanalyzer, topic: versions
    tuple val("${task.process}"), val('tabix'),            eval('tabix -h 2>&1 | grep Version | cut -f 2 -d " "')    , emit: versions_tabix           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: '' // MultiVCFAnalyzer has strict and input ordering and all are mandatory. Deactivating $args to prevent breakage of input
    def args2 = task.ext.args2 ?: ''

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
    """
    stub:
    """
    echo "" | gzip > fullAlignment.fasta.gz
    touch info.txt
    echo "" | gzip > snpAlignment.fasta.gz
    echo "" | gzip > snpAlignmentIncludingRefGenome.fasta.gz
    touch snpStatistics.tsv
    touch snpTable.tsv
    touch snpTableForSnpEff.tsv
    touch snpTableWithUncertaintyCalls.tsv
    touch structureGenotypes.tsv
    touch structureGenotypes_noMissingData-Columns.tsv
    touch MultiVCFAnalyzer.json
    """
}
