// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VEPRS_CACHEBUILDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE'
        : 'quay.io/biocontainers/YOUR-TOOL-HERE'}"

    input:
    tuple val(meta), path(annotation), path(fasta), path(fai)
    val species
    val assembly
    val release

    output:
    tuple val(meta), path("${prefix}"), emit: cache
    tuple val("${task.process}"), val('veprs'), eval("vep --version | sed 's/^vep //'"), topic: versions, emit: versions_veprs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Without an annotation, --auto fetches the species' whole release from Ensembl.
    def model = annotation.toString() ==~ /.*\.gff3?(\.gz)?$/ ? "--gff3 ${annotation}" : "--gtf ${annotation}"
    def build = annotation
        ? "--species ${species} --assembly ${assembly} ${model} --genome-fasta ${fasta} --output-dir ${prefix}/${species}/${release}_${assembly}"
        : "--auto ${species} ${assembly} --dir_cache ${prefix}"
    """
    vep-cache-builder \\
        --release ${release} \\
        ${build} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/${species}/${release}_${assembly}
    touch ${prefix}/${species}/${release}_${assembly}/info.json
    """
}
