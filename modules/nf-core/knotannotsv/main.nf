process KNOTANNOTSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/776546fb3850096f76abd8b63e0eac6b26807e629060e28c81e50ce056f3b8dd/data':
        'community.wave.seqera.io/library/perl-excel-writer-xlsx_perl-sort-key_perl-yaml-libyaml_git:a0a865d27eeceb13' }"

    input:
    tuple val(meta), path(annotsv_tsv)
    val(knot_out_xl)

    output:
    tuple val(meta), path("*.html"), emit: html, optional: true
    tuple val(meta), path("*.xlsm"), emit: xl, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def knot_prefix = task.ext.prefix ? "--outPrefix ${task.ext.prefix}" : "" // For knotAnnotSV, this a true prefix
    def knot_version = 'v1.1.5' // CHANGE when UPDATE
    def knot_script = knot_out_xl ? 'knotAnnotSV2XL.pl' : 'knotAnnotSV.pl'
    // TODO felix: Allow Excel output
    """
    git clone https://github.com/mobidic/knotAnnotSV.git --branch ${knot_version} --single-branch

    perl knotAnnotSV/${knot_script} \\
        ${args} \\
        --configFile knotAnnotSV/config_AnnotSV.yaml \\
        ${knot_prefix} \\
        --annotSVfile ${annotsv_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        knotAnnotSV: \$(echo ${knot_version})
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def knot_prefix = task.ext.prefix // For knotAnnotSV, this a true prefix
    def knot_version = 'v1.1.5' // CHANGE when UPDATE
    """
    echo $args

    touch ${knot_prefix}_${meta.id}.html
    touch ${knot_prefix}_${meta.id}.xlsm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        knotAnnotSV: \$(echo ${knot_version})
    END_VERSIONS
    """
}
