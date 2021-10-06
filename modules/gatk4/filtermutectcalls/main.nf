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

        conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
        } else {
            container "quay.io/biocontainers/gatk4:4.2.0.0--0"
        }

    input:
    tuple val(meta) , path(vcf) , path(tbi) , path(stats) , path(orientationbias) , path(segmentation) , path(contaminationfile) , val(contaminationest)
    path fasta
    path fastaidx
    path dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz")    , emit: filteredvcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: filteredtbi
    tuple val(meta), path("*.filteringStats.tsv") , emit: filteringstats
    path "versions.yml"                           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bias_command = ''
    def segmentation_command = ''
    def contamination_command = ''
    def bias_list = []
    def segmentation_list = []
    def contamination_list = []
    def useconfile = false

    useconfile = contaminationfile ? true : false
    orientationbias.each() {a -> bias_list.add(" --orientation-bias-artifact-priors " + a)}
    bias_command = bias_list.join(' ')
    segmentation.each() {a -> segmentation_list.add(" --tumor-segmentation " + a)}
    segmentation_command = segmentation_list.join(' ')
    if(useconfile){
        contaminationfile.each() {a -> contamination_list.add(" --contamination-table " + a)}
        contamination_command = contamination_list.join(' ')
    } else {
        contamination_command = contaminationest ? " --contamination-estimate ${contaminationest} " : ''
    }

    """
    gatk FilterMutectCalls \\
        -R $fasta \\
        -V $vcf \\
        $bias_command \\
        $segmentation_command \\
        $contamination_command \\
        -O ${prefix}.filtered.vcf.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
