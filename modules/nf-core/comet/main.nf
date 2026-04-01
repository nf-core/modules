/*
 * Perform spectrum identification using Comet.
 */
process COMET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/comet-ms:2026011--h9ee0642_0'
        : 'biocontainers/comet-ms:2026011--h9ee0642_0'}"

    input:
    tuple val(meta), path(mzml), path(fasta), path(comet_params)

    output:
    tuple val(meta), path("${prefix}*.comet.params"), emit: params
    tuple val(meta), path("${prefix}*.sqt"), emit: sqt, optional: true
    tuple val(meta), path("${prefix}*.txt"), emit: txt, optional: true
    tuple val(meta), path("${prefix}*.pep.xml"), emit: pepxml, optional: true
    tuple val(meta), path("${prefix}*.mzid"), emit: mzid, optional: true
    tuple val(meta), path("${prefix}*.pin"), emit: pin, optional: true
    tuple val("${task.process}"), val('comet'), eval("comet 2>&1 | head -2 | tail -1 | sed 's;.*\"\\(.*\\).*\";\\1;g'"), topic: versions, emit: versions_comet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    PARAMS_FILE="${prefix}.comet.params"

    if [ "${comet_params}" != "\${PARAMS_FILE}" ]; then
        # the param file has a different name than the required -> copy it to the required name
        cp ${comet_params} \${PARAMS_FILE}
    fi

    # adjust runtime parameters
    sed -i 's;database_name =.*;database_name = ${fasta};' \${PARAMS_FILE}
    sed -i "s;^num_threads.*;num_threads = ${task.cpus};" \${PARAMS_FILE}

    # run comet
    comet \\
        ${args} \\
        -P\${PARAMS_FILE} \\
        -N${prefix} \\
        ${mzml}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.comet.params
    touch ${prefix}.sqt
    touch ${prefix}.txt
    touch ${prefix}.pep.xml
    touch ${prefix}.mzid
    touch ${prefix}.pin
    """
}
