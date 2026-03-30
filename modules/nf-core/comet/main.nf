/*
 * Perform spectrum identification using Comet.
 */
process COMET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/comet-ms:2026.01.1--h9ee0642_0'
        : 'biocontainers/comet-ms:2026.01.1--h9ee0642_0'}"

    input:
    tuple val(meta), path(mzml), path(fasta), path(comet_params)

    output:
    tuple val(meta), path("*.comet.params"), emit: params
    tuple val(meta), path("*.sqt"), emit: sqt, optional: true
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    tuple val(meta), path("*.pep.xml"), emit: pepxml, optional: true
    tuple val(meta), path("*.mzid"), emit: mzid, optional: true
    tuple val(meta), path("*.pin"), emit: pin, optional: true
    tuple val("${task.process}"), val('comet'), eval("comet 2>&1 | head -2 | tail -1 | sed 's;.*\"\\(.*\\).*\";\\1;g'"), topic: versions, emit: versions_comet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def output_sqt = task.ext.output_sqt == null ? 0 : (task.ext.output_sqt ? 1 : 0)
    def output_txt = task.ext.output_txt == null ? 0 : (task.ext.output_txt ? 1 : 0)
    def output_pepxml = task.ext.output_pepxml == null ? 0 : (task.ext.output_pepxml ? 1 : 0)
    def output_mzidentml = task.ext.output_mzidentml == null ? 1 : (task.ext.output_mzidentml ? 1 : 0)
    def output_percolator = task.ext.output_percolator == null ? 0 : (task.ext.output_percolator ? 1 : 0)
    def num_output_lines = task.ext.num_output_lines == null ? 5 : task.ext.num_output_lines

    def comet_threads = 8
    if (!task.cpus) {
        log.info('[Comet] Available CPUs not known - defaulting to 8. Specify the number of CPUs to change this.')
    }
    else {
        comet_threads = task.cpus.intValue()
    }

    """
    BASE_FILENAME="${prefix}"
    PARAMS_FILE="${prefix}.comet.params"

    if [ "${comet_params}" != "\${PARAMS_FILE}" ]; then
        # the param file has a different name than the required -> copy it to the required name
        cp ${comet_params} \${PARAMS_FILE}
    fi

    # adjust runtime parameters
    sed -i 's;database_name =.*;database_name = ${fasta};' \${PARAMS_FILE}

    sed -i "s;^num_threads.*;num_threads = ${comet_threads};" \${PARAMS_FILE}

    sed -i "s;^output_sqtfile.*;output_sqtfile = ${output_sqt};" \${PARAMS_FILE}
    sed -i "s;^output_txtfile.*;output_txtfile = ${output_txt};" \${PARAMS_FILE}
    sed -i "s;^output_pepxmlfile.*;output_pepxmlfile =  ${output_pepxml};" \${PARAMS_FILE}
    sed -i "s;^output_mzidentmlfile.*;output_mzidentmlfile = ${output_mzidentml};" \${PARAMS_FILE}
    sed -i "s;^output_percolatorfile.*;output_percolatorfile = ${output_percolator};" \${PARAMS_FILE}

    sed -i "s;^num_output_lines.*;num_output_lines = ${num_output_lines};" \${PARAMS_FILE}

    # run comet
    comet \\
        ${args} \\
        -P\${PARAMS_FILE} \\
        -N\${BASE_FILENAME} \\
        ${mzml}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.comet.params

    touch ${prefix}.sqt
    touch ${prefix}.txt
    touch ${prefix}.pep.xml
    touch ${prefix}.mzid
    touch ${prefix}.pin
    """
}
