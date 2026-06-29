process DOTMATCH_PANEL_CHECK {
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
    tuple val(meta), path(panel)
    val k
    val metric

    output:
    tuple val(meta), path("panel_check"), emit: panel_check_dir
    tuple val(meta), path("panel_check/panel_summary.json"), emit: panel_summary
    tuple val(meta), path("panel_check/target_safety.tsv"), emit: target_safety
    tuple val(meta), path("panel_check/collision_pairs.tsv"), emit: collision_pairs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    dotmatch panel check ${panel} \\
      --k ${k} \\
      --metric ${metric} \\
      --out-dir panel_check \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
