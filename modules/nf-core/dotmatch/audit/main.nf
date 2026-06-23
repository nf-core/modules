process DOTMATCH_AUDIT {
    tag "$meta.id"
    label 'process_low'

    cpus   { task.ext.cpus   ?: 2 }
    memory { task.ext.memory ?: 4.GB }
    time   { task.ext.time   ?: 2.h  }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dotmatch:0.1.8--py314h118bc1c_0' :
        'quay.io/biocontainers/dotmatch:0.1.8--py314h118bc1c_0' }"

    input:
    tuple val(meta), path(targets)
    val k
    val audit_mode

    output:
    tuple val(meta), path("audit"), emit: audit_dir
    tuple val(meta), path("audit/audit_summary.json"), emit: summary_json
    tuple val(meta), path("audit/audit_summary.tsv"), emit: summary_tsv
    tuple val(meta), path("audit/target_safety.tsv"), emit: target_safety
    tuple val(meta), path("audit/collision_pairs.tsv"), emit: collision_pairs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    dotmatch audit \\
      --targets ${targets} \\
      --k ${k} \\
      --audit-mode ${audit_mode} \\
      --out-dir audit \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
