process DOTMATCH_ASSAY_RUN {
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
    tuple val(meta), path(assay_spec), path(assay_inputs)

    output:
    tuple val(meta), path("assay_report.html"), emit: assay_report
    tuple val(meta), path("assay_manifest.json"), emit: assay_manifest
    tuple val(meta), path("assay_manifest.summary.tsv"), emit: assay_manifest_summary
    tuple val(meta), path("sample_qc.tsv"), emit: sample_qc
    tuple val(meta), path("crispr_qc.html"), emit: crispr_qc_report
    tuple val(meta), path("crispr_qc.json"), emit: crispr_qc_json
    tuple val(meta), path("crispr_qc.summary.tsv"), emit: crispr_qc_summary
    tuple val(meta), path("counts.mageck.tsv"), emit: counts
    tuple val(meta), path("summary.json"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    for input_file in ${assay_inputs}; do
      input_base="\$(basename "\${input_file}")"
      cp -L "\${input_file}" ".\${input_base}.tmp"
      mv ".\${input_base}.tmp" "\${input_base}"
    done

    cp -L ${assay_spec} assay_spec.toml
    dotmatch assay run assay_spec.toml ${args}

    cp assay_out/assay_report.html assay_report.html
    cp assay_out/assay_manifest.json assay_manifest.json
    cp assay_out/assay_manifest.summary.tsv assay_manifest.summary.tsv
    cp assay_out/sample_qc.tsv sample_qc.tsv
    cp assay_out/crispr_qc.html crispr_qc.html
    cp assay_out/crispr_qc.json crispr_qc.json
    cp assay_out/crispr_qc.summary.tsv crispr_qc.summary.tsv
    cp assay_out/counts.mageck.tsv counts.mageck.tsv
    cp assay_out/summary.json summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
