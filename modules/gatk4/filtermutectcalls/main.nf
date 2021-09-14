// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta) , path(vcf) , path(tbi)

    path fasta
    path fastaidx
    path dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), emit: filteredvcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: filteredtbi
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk FilterMutectCalls \\
        -R $fasta \\
        -V $vcf \\
        -O ${prefix}.filtered.vcf.gz \\
        $options.args

    echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > ${software}.version.txt
    """
}
