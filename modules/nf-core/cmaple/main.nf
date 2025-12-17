process CMAPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmaple:1.1.0--h503566f_1':
        'biocontainers/cmaple:1.1.0--h503566f_1' }"

    input:
    tuple val(meta), path(aln), path(newick)

    output:
    tuple val(meta), path("*.treefile"), emit: treefile
    tuple val(meta), path("*.log")     , emit: log
    tuple val("${task.process}"), val("cmaple"), eval('cmaple --help | grep -m1 -oE "[0-9]+(\\.[0-9]+)+"'), topic: versions, emit: versions_cmaple

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def is_compressed    = aln.getExtension() == "gz"
    def aln_name         = is_compressed ? aln.getBaseName() : aln
    def uncompress_input = is_compressed ? "gzip -c -d ${aln} > ${aln_name}" : ''
    def tree_arg         = newick ? "-t ${newick}" : ""
    """
    $uncompress_input

    cmaple-aa \\
        $args \\
        -nt $task.cpus \\
        --prefix ${prefix} \\
        ${tree_arg} \\
        -aln $aln_name
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.treefile
    touch ${prefix}.log
    """
}
