// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MANTA_GERMLINE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::manta=1.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1"
    } else {
        container "quay.io/biocontainers/manta:1.6.0--h9ee0642_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path target_bed
    path target_bed_tbi

    output:
    tuple val(meta), path("*candidateSV.vcf.gz")             , emit: vcf
    tuple val(meta), path("*candidateSV.vcf.gz.tbi")         , emit: vcf_tbi
    tuple val(meta), path("*candidateSmallIndels.vcf.gz")    , emit: csi_vcf
    tuple val(meta), path("*candidateSmallIndels.vcf.gz.tbi"), emit: csi_vcf_tbi
    tuple val(meta), path("*tumorSV.vcf.gz")                 , emit: tumor_sv_vcf
    tuple val(meta), path("*tumorSV.vcf.gz.tbi")             , emit: tumor_sv_vcf_tbi
    path "versions.yml"                                      , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def options_manta = target_bed ? "--exome --callRegions $target_bed" : ""
    """
    configManta.py \
        --tumorBam $bam \
        --reference $fasta \
        $options_manta \
        --runDir manta

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
        ${prefix}.candidateSmallIndels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        ${prefix}.candidateSmallIndels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \
        ${prefix}.candidateSV.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
        ${prefix}.candidateSV.vcf.gz.tbi
    mv manta/results/variants/tumorSV.vcf.gz \
        ${prefix}.tumorSV.vcf.gz
    mv manta/results/variants/tumorSV.vcf.gz.tbi \
        ${prefix}.tumorSV.vcf.gz.tbi


    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( configManta.py --version )
    END_VERSIONS
    """
}
