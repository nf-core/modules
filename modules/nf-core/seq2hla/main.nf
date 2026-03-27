process SEQ2HLA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c8c9c11302a03e4f8fc9bbd0ce818e9f798e38f872c3ec099b8bb1c444185e0/data'
        : 'community.wave.seqera.io/library/seq2hla:2.3--26344ebe67dd0e8e'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*ClassI-class.bowtielog")            , emit: class1_bowtielog
    tuple val(meta), path("*ClassI-class.expression")           , emit: class1_expression
    tuple val(meta), path("*ClassI-class.HLAgenotype2digits")   , emit: class1_genotype_2d
    tuple val(meta), path("*ClassI-class.HLAgenotype4digits")   , emit: class1_genotype_4d
    tuple val(meta), path("*ClassI-nonclass.bowtielog")         , emit: class1_nonclass_bowtielog
    tuple val(meta), path("*ClassI-nonclass.expression")        , emit: class1_nonclass_expression
    tuple val(meta), path("*ClassI-nonclass.HLAgenotype2digits"), emit: class1_nonclass_genotype_2d
    tuple val(meta), path("*ClassI-nonclass.HLAgenotype4digits"), emit: class1_nonclass_genotype_4d
    tuple val(meta), path("*ClassII.bowtielog")                 , emit: class2_bowtielog
    tuple val(meta), path("*ClassII.expression")                , emit: class2_expression
    tuple val(meta), path("*ClassII.HLAgenotype2digits")        , emit: class2_genotype_2d
    tuple val(meta), path("*ClassII.HLAgenotype4digits")        , emit: class2_genotype_4d
    tuple val(meta), path("*.ambiguity")                        , emit: ambiguity, optional: true
    tuple val("${task.process}"), val('seq2hla'), eval("seq2HLA --version | sed 's/seq2HLA.py //'"), topic: versions, emit: versions_seq2hla

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

    seq2HLA \\
        -1 ${prefix}_1.fastq.gz \\
        -2 ${prefix}_2.fastq.gz \\
        -r ${prefix} \\
        -p ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-ClassI-class.HLAgenotype2digits
    touch ${prefix}-ClassII.HLAgenotype2digits
    touch ${prefix}-ClassI-class.HLAgenotype4digits
    touch ${prefix}-ClassII.HLAgenotype4digits
    touch ${prefix}-ClassI-class.bowtielog
    touch ${prefix}-ClassII.bowtielog
    touch ${prefix}-ClassI-class.expression
    touch ${prefix}-ClassII.expression
    touch ${prefix}-ClassI-nonclass.HLAgenotype2digits
    touch ${prefix}-ClassI-nonclass.HLAgenotype4digits
    touch ${prefix}-ClassI-nonclass.bowtielog
    touch ${prefix}-ClassI-nonclass.expression
    touch ${prefix}.ambiguity
    """
}
