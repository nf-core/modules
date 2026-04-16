
process HUMANN_RENORM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::humann=3.0.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/humann:3.0.0--pyh5e36f6f_1'
        : 'quay.io/biocontainers/humann:3.0.0--pyh5e36f6f_1'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_renorm.tsv.gz"), emit: renorm
    tuple val("${task.process}"), val('HUMAnN'), eval("humann --version 2>&1 | sed 's/humann v//'"), emit: versions_humann, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ $input == *.gz ]]; then
        gunzip -c $input > input.tsv
    else
        mv $input input.tsv
    fi

    humann_renorm_table \\
        --input input.tsv \\
        --output ${prefix}_renorm.tsv \\
        $args

    gzip -n ${prefix}_renorm.tsv

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "stub" | gzip > ${prefix}_renorm.tsv.gz

    """
}
