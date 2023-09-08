process KRAKENUNIQ_BUILD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::krakenuniq=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.2--pl5321h19e8d03_0':
        'biocontainers/krakenuniq:1.0.2--pl5321h19e8d03_0' }"

    input:
    tuple val(meta), path(custom_library_dir, stageAs: "library/*"), path(custom_taxonomy_dir, stageAs: "taxonomy"), path(custom_seqid2taxid)

    output:
    tuple val(meta), path("$prefix/"), emit: db
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    custom_db     = custom_library_dir ? "mkdir $prefix && mv library taxonomy $custom_seqid2taxid $prefix" : ""
    """
    $custom_db

    krakenuniq-build \\
        $args \\
        --threads ${task.cpus} \\
        --db ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """
}
