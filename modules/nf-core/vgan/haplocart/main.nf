process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:3.0.0--h9ee0642_0':
        'biocontainers/vgan:3.0.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    val(is_interleaved)

    output:
    tuple val(meta), path("*.output.txt")    , emit: txt
    tuple val(meta), path("*.posterior.txt") , emit: posterior
    tuple val("${task.process}"), val('vgan'), eval("vgan version 2>&1 | sed 's/vgan version //;s/ (Mela)//'"), topic: versions, emit: versions_vgan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_args = (meta.single_end || is_interleaved ) ? "-fq1 ${reads}" : "-fq1 ${reads[0]} -fq2 ${reads[1]}"
    def interleaved = is_interleaved ? "-i" : ""
    // on docker/singularity, need to copy the hcfiles as they lack some permissions
    def copy_hcfiles = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0 ? "cp -r /usr/local/share/vgan/hcfiles tmp_hcfiles" : ""
    def use_copied_hcfiles = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0 ? "--hcfiles tmp_hcfiles" : ""
    """
    ${copy_hcfiles}
    chmod 777 tmp_hcfiles
    chmod 777 tmp_hcfiles/*

    vgan haplocart \\
        $args \\
        -t $task.cpus \\
        $reads_args \\
        $interleaved \\
        -o ${prefix}.output.txt \\
        -pf ${prefix}.posterior.txt \\
        ${use_copied_hcfiles}
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}.output.txt
    touch ${prefix}.posterior.txt
    """
}
