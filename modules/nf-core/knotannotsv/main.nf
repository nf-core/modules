process KNOTANNOTSV {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/knotannotsv:1.1.5--hdfd78af_0'
        : 'biocontainers/knotannotsv:1.1.5--hdfd78af_0'}"

    input:
    tuple val(meta), path(annotsv_tsv), val(knot_out_xl)

    output:
    tuple val(meta), path("*.html"), emit: html, optional: true
    tuple val(meta), path("*.xlsm"), emit: xl, optional: true
    // CHANGE bellow when UPDATE
    tuple val("${task.process}"), val('knotannotsv'), val('1.1.5'), emit: versions_knotannotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Bellow for knotAnnotSV, this a true prefix
    prefix = task.ext.prefix ?: "${meta.id}"
    def knot_prefix = task.ext.prefix ? "--outPrefix ${task.ext.prefix}" : ""
    def knot_script = knot_out_xl ? 'knotAnnotSV2XL.pl' : 'knotAnnotSV.pl'
    def config_file = "\${CONDA_PREFIX:-/usr/local}/share/knotAnnotSV/config_AnnotSV.yaml"
    """
    ${knot_script} \\
        ${args} \\
        --configFile ${config_file} \\
        ${knot_prefix} \\
        --annotSVfile ${annotsv_tsv}
    """

    stub:
    def args = task.ext.args ?: ''
    // Bellow for knotAnnotSV, this a true prefix
    def knot_prefix = task.ext.prefix ? "${task.ext.prefix}_" : ""
    """
    echo ${args}

    touch ${knot_prefix}${meta.id}.html
    touch ${knot_prefix}${meta.id}.xlsm
    """
}
