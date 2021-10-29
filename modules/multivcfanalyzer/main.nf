// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIVCFANALYZER {
    tag '$bam'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }
    
    conda (params.enable_conda ? "bioconda::multivcfanalyzer=0.85.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multivcfanalyzer:0.85.2--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/multivcfanalyzer:0.85.2--hdfd78af_1"
    }

    input:
    path(vcfs)
    path(reference)
    val(allele_freq)
    val(min_qual)
    val(min_cov)
    val(min_freq_homo)
    val(min_freq_hetero)
    path snp_eff_optional
    path reference_gff_optional
    path pos_to_ignore_optional

    output:
    path "versions.yml", emit: versions
    path "fullAlignment.fasta", emit: fullalignment
    path "info.txt", emit: info
    path "snpAlignment.fasta", emit: snpalignment
    path "snpAlignmentIncludingRefGenome.fasta", emit: snpalignmentref
    path "snpStatistics.tsv", emit: snpstatistics
    path "snpTableForSnpeEff.tsv", emit: snptableforsnpeff
    path "snpTable.tsv", emit: snptable
    path "snpTableWithUncertaintyCalls.tsv", emit: snptablewithuncertaintycalls
    path "structureGenoypes_noMissingData-Columns.tsv", emit: structuregenotypesnomissing
    path "structureGenotypes.tsv", emit: structuregenotypes

    script:
    def vcf_list = vcfs.join(" ") : "$vcfs"
    def allele_freqs = !allele_freq ? log.error '[MULTIVCFANALYZER] required parameter of writing allele frequencies (vals T or F) not specified' : "$allele_freq"
    def min_quals = !min_qual ? log.error '[MULTIVCFANALYZER] required parameter of minimum genotyping quality not specified' : "$min_qual"
    def min_covs = !min_cov ? log.error '[MULTIVCFANALYZER] required parameter of minimum coverage for base call not specified' : "$min_cov"
    def min_freq_homos = !min_freq_homo ? log.error '[MULTIVCFANALYZER] required parameter of minimal allele frequency for homozygous call not specified' : "$min_freq_homo"
    def min_freq_heteros = !min_freq_hetero ? log.error '[MULTIVCFANALYZER] required parameter of minimal allele frequency for heterozygous call not specified' : "$min_freq_hetero"
    def snp_eff = snp_eff_optional ? $snp_eff_optional : "NA"
    def reference_gff = snp_eff_optional ? $snp_eff_optional : "NA"
    def pos_to_ignore = pos_to_ignore_optional ? $pos_to_ignore_optional : "NA"
    
    if (!task.memory) {
    log.info 'Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
    avail_mem = task.memory.giga
    }
    
    """ 

    multivcfanalyzer \\
        Xmx${avail_mem}g \\
        $snp_eff \\
        $reference \\
        $reference_gff \\
        ./ \\
        $allele_freqs \\
        $min_quals \\
        $min_covs \\
        $min_freq_homos \\
        $min_freq_heteros \\
        $pos_to_ignore \\
        $vcf_list

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( multivcfanalyzer --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """
}
