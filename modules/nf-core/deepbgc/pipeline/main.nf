process DEEPBGC_PIPELINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    input:
    tuple val(meta), path(genome)
    path db

    output:
    tuple val(meta), path("${prefix}/README.txt"), optional: true, emit: readme
    tuple val(meta), path("${prefix}/LOG.txt"), emit: log
    tuple val(meta), path("${prefix}/${prefix}.antismash.json"), optional: true, emit: json
    tuple val(meta), path("${prefix}/${prefix}.bgc.gbk"), optional: true, emit: bgc_gbk
    tuple val(meta), path("${prefix}/${prefix}.bgc.tsv"), optional: true, emit: bgc_tsv
    tuple val(meta), path("${prefix}/${prefix}.full.gbk"), optional: true, emit: full_gbk
    tuple val(meta), path("${prefix}/${prefix}.pfam.tsv"), optional: true, emit: pfam_tsv
    tuple val(meta), path("${prefix}/evaluation/${prefix}.bgc.png"), optional: true, emit: bgc_png
    tuple val(meta), path("${prefix}/evaluation/${prefix}.pr.png"), optional: true, emit: pr_png
    tuple val(meta), path("${prefix}/evaluation/${prefix}.roc.png"), optional: true, emit: roc_png
    tuple val(meta), path("${prefix}/evaluation/${prefix}.score.png"), optional: true, emit: score_png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export DEEPBGC_DOWNLOADS_DIR=${db}

    deepbgc \\
        pipeline \\
        ${args} \\
        ${genome}

    if [[ "${genome.baseName}/" != "${prefix}/" ]]; then
        mv "${genome.baseName}/" "${prefix}/"
    fi

    for i in \$(find -name '${genome.baseName}*' -type f); do
        mv \$i \${i/${genome.baseName}/${prefix}};
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(echo \$(deepbgc info 2>&1 /dev/null/ | grep 'version' | cut -d " " -f3) )
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/evaluation
    touch ${prefix}/README.txt
    touch ${prefix}/LOG.txt
    touch ${prefix}/${prefix}.antismash.json
    touch ${prefix}/${prefix}.bgc.gbk
    touch ${prefix}/${prefix}.bgc.tsv
    touch ${prefix}/${prefix}.full.gbk
    touch ${prefix}/${prefix}.pfam.tsv
    touch ${prefix}/evaluation/${prefix}.bgc.png
    touch ${prefix}/evaluation/${prefix}.pr.png
    touch ${prefix}/evaluation/${prefix}.roc.png
    touch ${prefix}/evaluation/${prefix}.score.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(echo \$(deepbgc info 2>&1 /dev/null/ | grep 'version' | cut -d " " -f3) )
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """
}
