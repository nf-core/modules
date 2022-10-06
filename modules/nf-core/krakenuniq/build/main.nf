
process KRAKENUNIQ_BUILD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::krakenuniq=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.0--pl5321h19e8d03_0':
        'quay.io/biocontainers/krakenuniq:1.0.0--pl5321h19e8d03_0' }"

    input:
    tuple path(custom_library_dir, stageAs: "${prefix}/lineage"), path(custom_taxonomy_dir, stageAs: "${prefix}/taxonomy"), path(custom_seqid2taxid, stageAs: "${prefix}/*"),
    path

    output:
    tuple val(meta), path("*"), emit: db
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def lib_dir = custom_library_dir ? "--library-dir ${custom_library_dir}" : ""
    def tax_dir = custom_taxonomy_dir ? "--taxonomy-dir ${custom_taxonomy_dir}" : ""
    """
    krakenuniq-build \\
        $args \\
        --threads ${task.cpus} \\
        ${lib_dir} \\
        ${tax_dir} \\
        --db ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """
}
