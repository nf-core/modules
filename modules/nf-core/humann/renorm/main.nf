include { getConda; getContainer; getHumannVersion } from '../utils'

process HUMANN_RENORM {
    tag "$meta.id"
    label 'process_low'

    conda { getConda(task.ext.version ?: getHumannVersion(task.process)) }
    container { getContainer(task.ext.version ?: getHumannVersion(task.process)) }

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_renorm.tsv.gz"), emit: renorm
    tuple val("${task.process}"), val('HUMAnN'), eval("humann --version 2>&1 | sed 's/humann v//'"), emit: versions_humann, topic: versions

    when:
    task.ext.when == null || task.ext.when

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
    echo "" | gzip > ${prefix}_renorm.tsv.gz

    """
}
