// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTGA_FASTGA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/685e31663606feff0225f6e4b79ce332c8516716a5fad8f2a411f0671f1c9126/data'
        : 'community.wave.seqera.io/library/fastga:1.5.20260729--8cb4c20b64de3563'}"

    input:
    tuple val(meta), path(fasta1), path(fasta2)
    val output_format

    output:
    tuple val(meta), path("*.paf"), emit: paf, optional: true
    tuple val(meta), path("*.psl"), emit: psl, optional: true
    tuple val(meta), path("*.1aln"), emit: onealn, optional: true
    tuple val(meta), path("*.1gdb"), emit: onegdb, optional: true
    tuple val(meta), path("*.gix"), emit: gix, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('fastga'), val("1.5.20260729"), topic: versions, emit: versions_fastga

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_input = fasta2 ? "${fasta1} ${fasta2}" : "${fasta1}"
    def output_arg = ""
    def fmt_arg = ""
    def valid_paf_formats = ["paf", "pafx", "pafm", "pafxS", "pafxs", "pafmS", "pafms"]
    if (output_format in valid_paf_formats) {
        fmt_arg = "-${output_format}"
        output_arg = "> ${prefix}.paf"
    }
    else if (output_format == "psl") {
        fmt_arg = "-psl"
        output_arg = "> ${prefix}.psl"
    }
    else if (output_format == "1aln") {
        output_arg = "-1:${prefix}.1aln"
    }
    else {
        error("Unsupported output format: ${output_format}")
    }
    """
    FastGA \\
        -T${task.cpus} \\
        -v \\
        ${args} \\
        ${fmt_arg} \\
        ${fasta_input} \\
        ${output_arg} \\
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf
    touch ${prefix}.psl
    touch ${prefix}.1aln
    """
}
