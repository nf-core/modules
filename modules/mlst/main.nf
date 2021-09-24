// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MLST {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::mlst=2.19.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mlst:2.19.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/mlst:2.19.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mlst \\
      --threads $task.cpus \\
      $fasta \\
      > ${fasta}.tsv

    cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
        mlst: \$( mlst --version | grep -o '[0-9.]*' )
    END_VERSIONS

    """

}
