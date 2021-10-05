// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_LEARNREADORIENTATIONMODEL {
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
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.artifact-prior.tar.gz"), emit: artifactprior
    path "versions.yml"                             , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputsList = []
    def inputsCommand = ''

    f1r2.each() {a -> inputsList.add(" -I " + a ) }
    inputsCommand = inputsList.join( ' ' )

    """
    gatk LearnReadOrientationModel \\
        ${inputsCommand} \\
        -O ${prefix}.artifact-prior.tar.gz \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
