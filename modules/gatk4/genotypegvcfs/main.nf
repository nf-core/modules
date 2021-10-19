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
    tuple val(meta), path(gvcf), path(gvcfIndex)
    path  fasta
    path  fastaIndex
    path  fastaDict
    path  dbsnp
    path  dbsnpIndex
    path  intervalsBed

    output:
    tuple val(meta), path("*.genotyped.vcf.gz"), emit: vcf
    path  "versions.yml"                       , emit: versions

    script:
    def prefix          = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def dbsnpOption     = dbsnp ? "-D ${dbsnp}" : ""
    def intervalsOption = intervalsBed ? "-L ${intervalsBed}" : ""
    """
    gatk \\
        GenotypeGVCFs \\
        $options.args \\
        $intervalsOption \\
        $dbsnpOption \\
        -R $fasta \\
        -V $gvcf \\
        -O ${prefix}.genotyped.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
