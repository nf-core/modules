process DEEPBGC_PIPELINE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::deepbgc=0.1.30" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deepbgc:0.1.30--pyhb7b1952_1':
        'quay.io/biocontainers/deepbgc:0.1.30--pyhb7b1952_1' }"

    input:
    tuple val(meta), path(genome)
    path(db)

    output:
    tuple val(meta), path("${genome.baseName}/README.txt")                                 ,   optional: true, emit: readme
    tuple val(meta), path("${genome.baseName}/LOG.txt")                                    ,   emit: log
    tuple val(meta), path("${genome.baseName}/${genome.baseName}.antismash.json")          ,   optional: true, emit: json
    tuple val(meta), path("${genome.baseName}/${genome.baseName}.bgc.gbk")                 ,   optional: true, emit: bgc_gbk
    tuple val(meta), path("${genome.baseName}/${genome.baseName}.bgc.tsv")                 ,   optional: true, emit: bgc_tsv
    tuple val(meta), path("${genome.baseName}/${genome.baseName}.full.gbk")                ,   optional: true, emit: full_gbk
    tuple val(meta), path("${genome.baseName}/${genome.baseName}.pfam.tsv")                ,   optional: true, emit: pfam_tsv
    tuple val(meta), path("${genome.baseName}/evaluation/${genome.baseName}.bgc.png")      ,   optional: true, emit: bgc_png
    tuple val(meta), path("${genome.baseName}/evaluation/${genome.baseName}.pr.png")       ,   optional: true, emit: pr_png
    tuple val(meta), path("${genome.baseName}/evaluation/${genome.baseName}.roc.png")      ,   optional: true, emit: roc_png
    tuple val(meta), path("${genome.baseName}/evaluation/${genome.baseName}.score.png")    ,   optional: true, emit: score_png
    path "versions.yml"                                                                    ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export DEEPBGC_DOWNLOADS_DIR=${db}

    deepbgc \\
        pipeline \\
        $args \\
        $genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(echo \$(deepbgc info 2>&1 /dev/null/ | grep 'version' | cut -d " " -f3) )
    END_VERSIONS
    """
}
