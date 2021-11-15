// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COOLER_CLOAD {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0"
    }

    input:
    tuple val(meta), path(pairs), path(index)
    val cool_bin
    path chromsizes

    output:
    tuple val(meta), val(cool_bin), path("*.cool"), emit: cool
    path "versions.yml"                           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def nproc  = options.args.contains('pairix') || options.args.contains('tabix')? "--nproc ${task.cpus}" : ''

    """
    cooler cload \\
        $options.args \\
        $nproc \\
        ${chromsizes}:${cool_bin} \\
        $pairs \\
        ${prefix}.${cool_bin}.cool

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
