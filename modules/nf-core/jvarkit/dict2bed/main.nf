process JVARKIT_DICT2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'quay.io/biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta), path(dict_files)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('jvarkit'), eval("jvarkit -v"), emit: versions_jvarkit, topic: versions


    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: '' /* give a chance to run a command like '| cut -f1,2,3 |sort | uniq' */
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir TMP

    jvarkit -Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP \\
        dict2bed \\
        ${args} ${dict_files} \\
        ${args2} > ${prefix}.bed

    rm -rf TMP
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.bed"
    """
}
