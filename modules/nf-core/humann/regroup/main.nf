include { getConda; getContainer; getHumannVersion } from '../utils'

process HUMANN_REGROUP {
    tag "$meta.id"
    label 'process_low'

    conda { getConda(task.ext.version ?: getHumannVersion(task.process)) }
    container { getContainer(task.ext.version ?: getHumannVersion(task.process)) }

    input:
    tuple val(meta), path(input)
    val groups
    path utility_db

    output:
    tuple val(meta), path("*_regroup.tsv.gz"), emit: regroup
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
    STATIC_CONFIG=`python -c "import humann; print(humann.__file__.replace('__init__.py', 'humann.cfg'))"`
    cat \$STATIC_CONFIG  | sed "s|utility_mapping = .*|utility_mapping = ${utility_db}|g" > humann.cfg
    export HUMANN_CONFIG=humann.cfg
    humann_config --print
    humann_regroup_table \\
        --input input.tsv \\
        --output ${prefix}_regroup.tsv \\
        --groups $groups \\
        $args

    gzip -n ${prefix}_regroup.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_regroup.tsv.gz
    """
}
