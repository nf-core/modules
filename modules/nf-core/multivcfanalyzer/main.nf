process MULTIVCFANALYZER {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::multivcfanalyzer=0.85.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multivcfanalyzer:0.85.2--hdfd78af_1':
        'biocontainers/multivcfanalyzer:0.85.2--hdfd78af_1' }"

    input:
    path vcfs
    path fasta
    path snpeff_results
    path gff
    val allele_freqs
    val genotype_quality
    val coverage
    val homozygous_freq
    val heterozygous_freq
    path gff_exclude


    output:
    path('fullAlignment.fasta.gz')                       , emit: full_alignment
    path('info.txt')                                     , emit: info_txt
    path('snpAlignment.fasta.gz')                        , emit: snp_alignment
    path('snpAlignmentIncludingRefGenome.fasta.gz')      , emit: snp_genome_alignment
    path('snpStatistics.tsv')                            , emit: snpstatistics
    path('snpTable.tsv')                                 , emit: snptable
    path('snpTableForSnpEff.tsv')                        , emit: snptable_snpeff
    path('snpTableWithUncertaintyCalls.tsv')             , emit: snptable_uncertainty
    path('structureGenotypes.tsv')                       , emit: structure_genotypes
    path('structureGenotypes_noMissingData-Columns.tsv') , emit: structure_genotypes_nomissing
    path('MultiVCFAnalyzer.json')                        , emit: json
    path "versions.yml"                                  , emit: versions

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
        ${vcfs.join(" ")}

    gzip \\
        $args2 \\
        fullAlignment.fasta snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multivcfanalyzer: \$(echo \$(multivcfanalyzer --help | head -n 1) | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
