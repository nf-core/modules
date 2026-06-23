process DOTMATCH_CRISPR_COUNT {
    tag "$meta.id"
    label 'process_medium'

    cpus   { task.ext.cpus   ?: 4 }
    memory { task.ext.memory ?: 8.GB }
    time   { task.ext.time   ?: 4.h  }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dotmatch:0.1.8--py314h118bc1c_0' :
        'quay.io/biocontainers/dotmatch:0.1.8--py314h118bc1c_0' }"

    input:
    tuple val(meta), path(reads), path(library)
    val guide_start
    val guide_length
    val k
    val metric

    output:
    tuple val(meta), path("*.counts.mageck.tsv"), emit: counts
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.sample_qc.tsv"), emit: sample_qc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def threads = task.cpus ?: 1
    """
    printf 'sample_id\\tfastq\\n' > samples.tsv
    printf '%s\\t%s\\n' '${meta.id}' '${reads}' >> samples.tsv

    dotmatch crispr-count \\
      --library ${library} \\
      --samples samples.tsv \\
      --guide-start ${guide_start} \\
      --guide-length ${guide_length} \\
      --k ${k} \\
      --metric ${metric} \\
      --ambiguity-policy radius \\
      --threads ${threads} \\
      --out ${prefix}.counts.mageck.tsv \\
      --summary ${prefix}.summary.json \\
      --sample-qc ${prefix}.sample_qc.tsv \\
      --ambiguous discard \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
