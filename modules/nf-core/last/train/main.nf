process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1595--h43eeafb_0' :
        'biocontainers/last:1595--h43eeafb_0' }"

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.train"), emit: param_file
    tuple val(meta), path("*.tsv")  , emit: multiqc
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)

    last-train \\
        $args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        > ${prefix}.train

    echo "id\tsubstitution_percent_identity\tlast -t\tlast -a\tlast -A\tlast -b\tlast -B\tlast -S"       > ${prefix}.train.tsv
    printf "\$(basename ${prefix}.train .target.train)\t"                                               >> ${prefix}.train.tsv
    grep 'substitution percent identity' ${prefix}.train | tail -n 1 | awk '{print \$5}' | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -t' ${prefix}.train | tail -n 1 | awk '{print \$2}'   | sed -e 's/-t//'   | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -a' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -A' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -b' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -B' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\n' '\t' >> ${prefix}.train.tsv
    grep 'last -S' ${prefix}.train | tail -n 1 | awk '{print \$3}'                                      >> ${prefix}.train.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)
    touch ${prefix}.train
    touch ${prefix}.train.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version | sed 's/lastdb //')
    END_VERSIONS
    """
}
