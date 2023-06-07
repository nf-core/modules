process METABULI_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::metabuli=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.0--pl5321hf1761c0_0':
        'biocontainers/metabuli:1.0.0--pl5321hf1761c0_0' }"

    input:
    tuple val(meta), path(fastas)
    path(db)

    output:
    tuple val(meta), path("*/*_classifications.tsv"), emit: classification
    tuple val(meta), path("*/*_report.tsv"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = meta.single_end ? fastas.getName().endsWith(".gz") : fastas[0].getName().endsWith(".gz")
    def input = meta.single_end ? "--seq-mode 1 ${fastas}" : "${fastas[0]} ${fastas[1]}"
    if (is_compressed && meta.single_end) {
      input = "--seq-mode 1 ${fastas.baseName}"
    } else if (is_compressed) {
      input =  "${fastas[0].baseName} ${fastas[1].baseName}"
    }
    
    """
    if [ "$is_compressed" == "true" ]; then
    gzip -d *.gz
    fi

    metabuli \\
        classify \\
        $args \\
        --threads $task.cpus \\
        ${input} \\
        ${db} \\
        ${prefix}_out \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli | grep Version | sed 's/^metabuli Version: //';))
    END_VERSIONS
    """
}
