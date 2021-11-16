// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_GENOTYPEGVCFS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)
    path  fasta
    path  fasta_index
    path  fasta_dict
    path  dbsnp
    path  dbsnp_index
    path  intervals_bed

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    script:
    def prefix           = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def dbsnp_options    = dbsnp ? "-D ${dbsnp}" : ""
    def interval_options = intervals_bed ? "-L ${intervals_bed}" : ""
    def gvcf_options     = gvcf.name.endsWith(".vcf") || gvcf.name.endsWith(".vcf.gz") ? "$gvcf" : "gendb://$gvcf"
    """
    gatk \\
        GenotypeGVCFs \\
        $options.args \\
        $interval_options \\
        $dbsnp_options \\
        -R $fasta \\
        -V $gvcf_options \\
        -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
