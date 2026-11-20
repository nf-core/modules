// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VEPRS_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE'
        : 'quay.io/biocontainers/YOUR-TOOL-HERE'}"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val genome
    val species
    val cache_version
    tuple val(meta2), path(cache)
    tuple val(meta3), path(fasta), path(fai)
    path extra_files

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.tab.gz"), emit: tab, optional: true
    tuple val(meta), path("${prefix}.json.gz"), emit: json, optional: true
    tuple val(meta), val("${task.process}"), val('veprs'), path("*.html"), topic: multiqc_files, emit: report, optional: true
    tuple val("${task.process}"), val('veprs'), eval("vep --version | sed 's/^vep //'"), topic: versions, emit: versions_veprs
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | sed -n 's/^Version: *//p'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def create_index = file_extension == "vcf" ? "tabix ${args2} ${prefix}.${file_extension}.gz" : ""
    """
    vep \\
        -i ${vcf} \\
        -o ${prefix}.${file_extension}.gz \\
        ${args} \\
        ${compress_cmd} \\
        ${reference} \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache \${PWD}/${cache} \\
        --fork ${task.cpus}

    ${create_index}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}.${file_extension}.gz_summary.html
    """
}
