process NGSCHECKMATE_VAFNCM {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' :
        'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' }"

    input:
    tuple val(meta), path(vafs)

    output:
    tuple val(meta), path("*.pdf")             , emit: pdf, optional: true
    tuple val(meta), path("*_corr_matrix.txt") , emit: corr_matrix
    tuple val(meta), path("*_all.txt")         , emit: all
    tuple val(meta), path("*_matched.txt")     , emit: matched
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # tool has a bug where it misses the final file, so add a dummy one.
    cp ${vafs[0]} zzzzzz.vaf
    vaf_ncm.py -I . -O . -N ${prefix} $args

    # remove the existance of the dummy file
    rm zzzzzz.vaf
    sed -i.bak "/zzzzzz/d" ${prefix}_all.txt

    # generate a file with all the samples that do match, for consistency with the bam mode (ngscheckmate/ncm)
    sed "/unmatched/d" ${prefix}_all.txt > ${prefix}_matched.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}
