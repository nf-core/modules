process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/05ff1dba1dcf86a28dd9d911ffcbc925175b370e77b5ea699ae421b8481556bc/data'
        : 'community.wave.seqera.io/library/last_gzip:40e65bc1d7d8624e'}"

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.train"), emit: param_file
    tuple val(meta), path("*.tsv")  , emit: multiqc
   // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

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

    echo "id\tsubstitution_percent_identity\tlast -t\tlast -a\tlast -A\tlast -b\tlast -B\tlast -S"         > ${prefix}.train.tsv
    printf "\$(basename ${prefix}.train .target.train)\t"                                                 >> ${prefix}.train.tsv
    # Do not take the last 'substitution percent identity' value; it is calculated from a matrix rounded to integers
    grep 'substitution percent identity' ${prefix}.train |
        tail -n 2 | head -n 1 | awk '{print \$5}'                                        | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -t' ${prefix}.train | tail -n 1 | awk '{print \$2}'   | sed -e 's/-t//'   | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -a' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -A' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -b' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -B' ${prefix}.train | tail -n 1 | awk '{print \$3}'                       | tr '\\n' '\\t' >> ${prefix}.train.tsv
    grep 'last -S' ${prefix}.train | tail -n 1 | awk '{print \$3}'                                        >> ${prefix}.train.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)
    touch ${prefix}.train
    touch ${prefix}.train.tsv
    """
}
