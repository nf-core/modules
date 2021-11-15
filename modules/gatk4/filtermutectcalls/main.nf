// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(orientationbias), path(segmentation), path(contaminationfile), val(contaminationest)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz")            , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")        , emit: tbi
    tuple val(meta), path("*.filteringStats.tsv"), emit: stats
    path "versions.yml"                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def orientationbias_options = ''
    if (orientationbias) {
        orientationbias_options = '--orientation-bias-artifact-priors ' + orientationbias.join(' --orientation-bias-artifact-priors ')
    }

    def segmentation_options = ''
    if (segmentation) {
        segmentation_options = '--tumor-segmentation ' + segmentation.join(' --tumor-segmentation ')
    }

    def contamination_options = contaminationest ? " --contamination-estimate ${contaminationest} " : ''
    if (contaminationfile) {
        contamination_options = '--contamination-table ' + contaminationfile.join(' --contamination-table ')
    }
    """
    gatk FilterMutectCalls \\
        -R $fasta \\
        -V $vcf \\
        $orientationbias_options \\
        $segmentation_options \\
        $contamination_options \\
        -O ${prefix}.vcf.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
