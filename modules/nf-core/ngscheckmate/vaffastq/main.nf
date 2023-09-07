process NGSCHECKMATE_VAFFASTQ {
    label 'process_single'

    conda "bioconda::ngscheckmate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' :
        'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' }"

    input:
    path(vafs)

    output:
    path "*.pdf"            , emit: pdf, optional: true
    path "*_corr_matrix.txt", emit: corr_matrix
    path "*_matched.txt"    , emit: matched
    path "*_all.txt"        , emit: all
    path "*.vcf"            , emit: vcf, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "NGSCheckMate"

    """
    vaf_ncm.py -I . -O . -N ${prefix} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${prefix}_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
