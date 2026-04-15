include { getConda; getContainer; getExt; getProcessName } from '../utils'
def getProcessNamePrefix(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0]
}
process HUMANNREGROUP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? { getConda(getProcessNamePrefix(task.process)) } : null)
    //if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container)	{
    //    container { "docker://" + getContainer(getProcessNamePrefix(task.process)) }
    //} else {
    container { getContainer(getProcessNamePrefix(task.process)) }
    //}
    input:
    tuple val(meta), path(input)
    val groups
    path utility_db

    output:
    tuple val(meta), path("*_regroup.tsv.gz"), emit: regroup
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "stub" | gzip >  ${prefix}_regroup.tsv.gz
    """
}
