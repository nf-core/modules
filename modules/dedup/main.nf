process DEDUP {
    tag "$meta.id"
    label 'process_low'

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
    def args = task.ext.args ?: ''
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    dedup \\
        -Xmx${task.memory.toGiga()}g  \\
        -i $bam \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(dedup --version 2>&1) | tail -n 1 | sed 's/.* v//')

    END_VERSIONS
    """
}
