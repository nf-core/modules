// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process RAVEN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    conda (params.enable_conda ? "bioconda::raven-assembler=1.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://quay.io/biocontainers/raven-assembler:1.6.1--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/raven-assembler:1.6.1--h2e03b76_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa")     , emit: gfa
    path "versions.yml"                , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    raven \\
        -t $task.cpus \\
        --graphical-fragment-assembly ${prefix}.gfa \\
        $options.args \\
        $reads | \\
        gzip -c > ${prefix}.fasta.gz
    
    gzip -c ${prefix}.gfa > ${prefix}.gfa.gz
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( raven --version )
    END_VERSIONS
    """
}
