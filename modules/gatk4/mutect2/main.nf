// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'
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
    tuple val(meta) , path(input) , path(input_index) , val(which_norm)
    val  run_single
    val  run_pon
    val  run_mito
    val  interval_label
    path fasta
    path fai
    path dict
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def panels_command = ''
    def normals_command = ''

    def inputs_command = '-I ' + input.join( ' -I ')

    if(run_pon) {
        panels_command = ''
        normals_command = ''

    } else if(run_single) {
        panels_command = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals"
        normals_command = ''

    } else if(run_mito){
        panels_command = "-L ${interval_label} --mitochondria-mode"
        normals_command = ''

    } else {
        panels_command = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals --f1r2-tar-gz ${prefix}.f1r2.tar.gz"
        normals_command = '-normal ' + which_norm.join( ' -normal ')
    }

    """
    gatk Mutect2 \\
        -R ${fasta} \\
        ${inputs_command} \\
        ${normals_command} \\
        ${panels_command} \\
        -O ${prefix}.vcf.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
