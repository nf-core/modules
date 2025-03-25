process DEDUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dedup:0.12.8--hdfd78af_1' :
        'biocontainers/dedup:0.12.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_rmdup.bam"), emit: bam     // _rmdup is hardcoded output from dedup
    tuple val(meta), path("*.json")     , emit: json
    tuple val(meta), path("*.hist")     , emit: hist
    tuple val(meta), path("*log")       , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    dedup \\
        -Xmx${task.memory.toGiga()}g  \\
        -i $bam \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dedup: \$( echo \$(dedup --version 2>&1) | tail -n 1 | sed 's/.* v//')

    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paired_end.dedup.json
    touch ${prefix}.paired_end.hist
    touch ${prefix}.paired_end.log
    touch ${prefix}.paired_end_rmdup.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dedup: \$( echo \$(dedup --version 2>&1) | tail -n 1 | sed 's/.* v//')

    END_VERSIONS
    """
}
