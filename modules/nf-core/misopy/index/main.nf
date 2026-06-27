process MISOPY_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27_0':
        'quay.io/biocontainers/misopy:0.5.4--py27_0' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("index"), emit: miso_index
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('misopy'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'misopy\').version)"'), topic: versions, emit: versions_misopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -v "^#" ${gff3} | \\
    awk 'BEGIN{OFS="\\t"} NF==9 && \$6~/^[0-9.]+\$/ || \$6=="." {print}' > filtered.gff3

    index_gff \\
        --index filtered.gff3 \\
        "index"
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p "index"
    touch "index/${prefix}.shelve"
    """
}
