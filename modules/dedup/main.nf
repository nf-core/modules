// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEDUP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::dedup=0.12.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/dedup:0.12.8--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/dedup:0.12.8--hdfd78af_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_rmdup.bam"), emit: bam     // _rmdup is hardcoded output from dedup
    tuple val(meta), path("*.json")     , emit: json
    tuple val(meta), path("*.hist")     , emit: hist
    tuple val(meta), path("*log")       , emit: log
    path "versions.yml"                 , emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    dedup \\
        -Xmx${task.memory.toGiga()}g  \\
        -i $bam \\
        -o . \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(dedup --version 2>&1) | tail -n 1 | sed 's/.* v//')

    END_VERSIONS
    """
}
