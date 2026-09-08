process SATIVAEPANG_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.9.3.4--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.9.3.4--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(alignment), path(taxonomy), val(taxcode)
    path(reftree)
    path(refmodel)

    output:
    tuple val(meta), path("*.refjson"), emit: refjson
    tuple val(meta), path("*.model")  , emit: model
    tuple val("${task.process}"), val('sativaepang'), eval("sed -n 's#.*share/sativa-epang-\\([0-9.]*\\)-.*#\\1#p' \$(command -v sativa-epang)"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -reftree/-refmodel hand off a tree built elsewhere (e.g. RAxML-NG): EPA-ng then
    // numbers the branches and the taxonomy map/node heights are computed in pure Python,
    // so no RAxML runs at all. Without them this stage falls back to sativa-epang's own
    // constrained RAxML search.
    def reftree_arg  = reftree  ? "-reftree ${reftree}"   : ''
    def refmodel_arg = refmodel ? "-refmodel ${refmodel}" : ''
    """
    sativa-epang \\
        -s ${alignment} \\
        -t ${taxonomy} \\
        -x ${taxcode} \\
        ${reftree_arg} \\
        ${refmodel_arg} \\
        -n ${prefix} \\
        -o . \\
        -T ${task.cpus} \\
        -stage reference \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.refjson
    touch ${prefix}.model
    """
}
