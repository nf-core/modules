process KNOTANNOTSV {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/knotannotsv:1.1.5--hdfd78af_0'
        : 'quay.io/biocontainers/knotannotsv:1.1.5--hdfd78af_0'}"

    input:
    tuple val(meta), path(annotsv_tsv), val(knot_out_xl)

    output:
    tuple val(meta), path("*.{html,xlsm}"), emit: output_file
    // CHANGE bellow when UPDATE
    tuple val("${task.process}"), val('knotannotsv'), val('1.1.5'), emit: versions_knotannotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def knot_script = knot_out_xl ? 'knotAnnotSV2XL.pl' : 'knotAnnotSV.pl'
    def config_file = "\${CONDA_PREFIX:-/usr/local}/share/knotAnnotSV/config_AnnotSV.yaml"
    """
    # This tool always name outFile based on inFile (eg: sample.tsv -> sample.html)
    # So rename inFile on-the-fly to re-enable control of outName by 'ext.prefix'
    mv ${annotsv_tsv} ${prefix}.tsv

    ${knot_script} \\
        ${args} \\
        --configFile ${config_file} \\
        --annotSVfile ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_file = knot_out_xl ? "${prefix}.xlsm" : "${prefix}.html"
    """
    echo ${args}

    touch ${out_file}
    """
}
