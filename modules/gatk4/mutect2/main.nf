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

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta) , path(bam) , path(bai) , val(which_norm)
    val run_single
    val run_pon
    path fasta
    path fastaidx
    path dict
    path germline_resource
    path germline_resource_idx
    path panel_of_normals
    path panel_of_normals_idx

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputsList = []
    def normalsList = []
    def inputsCommand = ''
    def panelsCommand = ''
    def normalsCommand = ''

    bam.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ')

    if(run_pon) {
        panelsCommand = ''
        normalsCommand = ''

    } else if(run_single) {
        panelsCommand = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals"
        normalsCommand = ''

    } else {
        panelsCommand = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals --f1r2-tar-gz ${prefix}.f1r2.tar.gz"
        which_norm.each() {a -> normalsList.add(" -normal " + a ) }
        normalsCommand = normalsList.join( ' ')
    }

    """
    gatk Mutect2 \\
        -R ${fasta} \\
        ${inputsCommand} \\
        ${normalsCommand} \\
        ${panelsCommand} \\
        -O ${prefix}.vcf.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        - ${getSoftwareName(task.process)}: \$(gatk --version 2>&1 | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
