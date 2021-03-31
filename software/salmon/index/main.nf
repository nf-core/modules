// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', publish_id:'') }

    conda (params.enable_conda ? "bioconda::salmon=1.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.4.0--hf69c8f4_0"
    } else {
        container "quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"
    }

    input:
    path transcript_fasta

    output:
    path "salmon"       , emit: index
    path "*.version.txt", emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    salmon \\
        index \\
        --threads $task.cpus \\
        -t $transcript_fasta \\
        $options.args \\
        -i salmon
    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}
