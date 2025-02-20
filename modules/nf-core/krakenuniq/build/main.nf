process KRAKENUNIQ_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2'
        : 'biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2'}"

    input:
    tuple val(meta), path(custom_library_dir, stageAs: "library/*"), path(custom_taxonomy_dir, stageAs: "taxonomy"), path(custom_seqid2taxid)
    val keep_intermediate

    output:
    tuple val(meta), path("${prefix}/"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    custom_db = custom_library_dir ? "mkdir ${prefix} && mv library taxonomy ${custom_seqid2taxid} ${prefix}" : ""
    run_cleanup = keep_intermediate ? "" : "find -L ${prefix} -type f -not -name \"*.kdb\" -type f -not -name \"*idx\" -not -name \"taxDB\" -delete"

    """
    ${custom_db}

    krakenuniq-build \\
        ${args} \\
        --threads ${task.cpus} \\
        --db ${prefix}

    ${run_cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    run_cleanup = keep_intermediate ? "" : "find -L ${prefix} -type f -not -name \"*.kdb\" -type f -not -name \"*idx\" -not -name \"taxDB\" -delete"
    """
    mkdir ${prefix}/
    touch ${prefix}/database-build.log
    touch ${prefix}/database.idx
    touch ${prefix}/database.jdb
    touch ${prefix}/database.kdb
    touch ${prefix}/database.kdb.counts
    touch ${prefix}/database.kraken.tsv
    touch ${prefix}/database.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakenuniq: \$(echo \$(krakenuniq --version 2>&1) | sed 's/^.*KrakenUniq version //; s/ .*\$//')
    END_VERSIONS
    """
}
