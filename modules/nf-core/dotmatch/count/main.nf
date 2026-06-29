process DOTMATCH_COUNT {
    tag "$meta.id"
    label 'process_medium'

    cpus   { task.ext.cpus   ?: 1 }
    memory { task.ext.memory ?: 8.GB }
    time   { task.ext.time   ?: 4.h  }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dotmatch:0.1.8--py314h118bc1c_0' :
        'quay.io/biocontainers/dotmatch:0.1.8--py314h118bc1c_0' }"

    input:
    tuple val(meta), path(reads), path(targets)
    val target_start
    val target_length
    val k
    val metric

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.sample_qc.tsv"), emit: sample_qc
    tuple val(meta), path("*.assignments.tsv"), emit: assignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def sampleLabel = meta.sample_label ?: meta.id
    def threads = task.cpus ?: 1
    """
    dotmatch count \\
      --targets ${targets} \\
      --reads ${reads} \\
      --sample-label ${sampleLabel} \\
      --target-start ${target_start} \\
      --target-length ${target_length} \\
      --k ${k} \\
      --metric ${metric} \\
      --ambiguity-policy radius \\
      --threads ${threads} \\
      --out ${prefix}.counts.tsv \\
      --summary ${prefix}.summary.json \\
      --sample-qc ${prefix}.sample_qc.tsv \\
      --assignments ${prefix}.assignments.tsv \\
      --ambiguous discard \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
