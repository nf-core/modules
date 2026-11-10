// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTVEP_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3dda8b0c7fd437d6a0e82050ea78debf20dfbcb566dd8b63dbfecef51ae5979/data'
        : 'community.wave.seqera.io/library/fastvep_htslib:1cb9b1d86ed7c379'}"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    tuple val(meta2), path(gff3)
    tuple val(meta3), path(fasta), path(fai)
    path(sa_files, stageAs: 'fastvep_sa/*')
    tuple val(meta4), path(cache_dir)
    val output_format

    output:
    tuple val(meta), path("${prefix}.vcf.gz")           , optional: true, emit: vcf
    tuple val(meta), path("${prefix}.tab.gz")           , optional: true, emit: tab
    tuple val(meta), path("${prefix}.json.gz")          , optional: true, emit: json
    tuple val(meta2), path("${gff3}.fastvep.cache")     , optional: true, emit: transcript_cache
    tuple val("${task.process}"), val('fastvep'), eval("fastvep --version | sed 's/^fastvep //'"), topic: versions, emit: versions_fastvep
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_annotated"
    if (!(output_format in ['vcf', 'tab', 'json'])) {
        error "FASTVEP_ANNOTATE: invalid output_format '${output_format}'. Must be one of: vcf, tab, json."
    }
    def fasta_arg = fasta ? "--fasta \"${fasta}\"" : ''
    def sa_arg = sa_files ? "--sa-dir fastvep_sa" : ''
    def cache_dir_arg = cache_dir ? "--cache-dir \"${cache_dir}\"" : ''
    """
    fastvep \\
        annotate \\
        --input "${vcf}" \\
        --output "${prefix}.${output_format}" \\
        --output-format ${output_format} \\
        --gff3 "${gff3}" \\
        ${fasta_arg} \\
        ${sa_arg} \\
        ${cache_dir_arg} \\
        --no-progress \\
        ${args}

    bgzip \\
        --threads ${task.cpus} \\
        ${args2} \\
        "${prefix}.${output_format}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_annotated"
    if (!(output_format in ['vcf', 'tab', 'json'])) {
        error "FASTVEP_ANNOTATE: invalid output_format '${output_format}'. Must be one of: vcf, tab, json."
    }
    """
    echo "" | bgzip > "${prefix}.${output_format}.gz"
    touch "${gff3}.fastvep.cache"
    """
}
