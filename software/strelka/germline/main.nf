// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
    } else {
        container "quay.io/biocontainers/strelka:2.9.10--0"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path target_bed

    output:
    tuple val(meta), path("*_variants.vcf.gz"), emit: vcf
    tuple val(meta), path("*_variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*_genome.vcf.gz"), emit: genome_vcf
    tuple val(meta), path("*_genome.vcf.gz.tbi"), emit: genome_vcf_tbi
    path "*.version.txt", emit: version 

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "strelka_${meta.id}${ioptions.suffix}" : "strelka_${meta.id}"
    
    options_strelka = params.target_bed ? "--exome --callRegions ${target_bed}" : ""
    """
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options_strelka} \
        ${options.args} \
        --runDir strelka
    python strelka/runWorkflow.py -m local -j ${task.cpus}
    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}_genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}_genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}_variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}_variants.vcf.gz.tbi
    echo configureStrelkaGermlineWorkflow.py --version &> ${software}.version.txt #2>&1
    """
}
