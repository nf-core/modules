process CNVPYTOR_CALLCNVS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.2.1--pyhdfd78af_0':
        'biocontainers/cnvpytor:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pytor)
    val bin_sizes

    output:
    tuple val(meta), path("${pytor.baseName}.pytor")	, emit: pytor
    tuple val("${task.process}"), val('cnvpytor'), eval("cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p'"), emit: versions_cnvpytor, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bins_cmd = bin_sizes ? "-call $bin_sizes" : '-call 1000'
    """
    cnvpytor \\
        -root ${pytor} \\
        ${bins_cmd}
    """

    stub:
    """
    touch ${pytor.baseName}.pytor
    """
}
